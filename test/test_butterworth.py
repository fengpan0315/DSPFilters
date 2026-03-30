import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import signal
from scipy.signal import sos2zpk
import re
import sys
import argparse
from pathlib import Path

class DspFiltersReader:
    def __init__(self, filename):
        self.filename = filename
        self.filter_type = ""
        self.filter_info = {}
        self.biquads = []
        self.read_file()

    def read_file(self):
        with open(self.filename, 'r') as f:
            lines = f.readlines()

        for line in lines:
            if not self.filter_type:
                if 'Lowpass' in line:
                    self.filter_type = 'lowpass'
                elif 'Highpass' in line:
                    self.filter_type = 'highpass'
                elif 'Bandpass' in line:
                    self.filter_type = 'bandpass'
                elif 'Bandstop' in line:
                    self.filter_type = 'bandstop'

            if 'Order:' in line:
                self.filter_info['order'] = int(line.split(':')[1].strip())
            elif 'Sample Rate:' in line:
                self.filter_info['sample_rate'] = float(line.split(':')[1].strip().split()[0])
            elif 'Cutoff Frequency:' in line:
                self.filter_info['cutoff'] = float(line.split(':')[1].strip().split()[0])
            elif 'Center Frequency:' in line:
                self.filter_info['center'] = float(line.split(':')[1].strip().split()[0])
            elif 'Bandwidth:' in line:
                self.filter_info['bandwidth'] = float(line.split(':')[1].strip().split()[0])

        pattern = r'Stage (\d+):\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)'

        for line in lines:
            match = re.search(pattern, line)
            if match:
                self.biquads.append({
                    'a0': float(match.group(2)),
                    'a1': float(match.group(3)),
                    'a2': float(match.group(4)),
                    'b0': float(match.group(5)),
                    'b1': float(match.group(6)),
                    'b2': float(match.group(7))
                })

    def get_sos(self):
        sos = []
        for bq in self.biquads:
            sos.append([bq['b0'], bq['b1'], bq['b2'],
                       bq['a0'], bq['a1'], bq['a2']])
        return np.array(sos)


class SciPyButterworth:
    def __init__(self, filter_type, order, sample_rate, **kwargs):
        self.filter_type = filter_type
        self.order = order
        self.sample_rate = sample_rate
        self.sos = self.design_filter(**kwargs)

    def design_filter(self, **kwargs):
        if self.filter_type in ['lowpass', 'highpass']:
            cutoff = kwargs.get('cutoff', 1000)
            nyquist = self.sample_rate / 2
            Wn = cutoff / nyquist
            return signal.butter(self.order, Wn, btype=self.filter_type,
                                output='sos')

        elif self.filter_type in ['bandpass', 'bandstop']:
            center = kwargs.get('center', 1000)
            bandwidth = kwargs.get('bandwidth', 100)
            nyquist = self.sample_rate / 2
            Wn = [center - bandwidth/2, center + bandwidth/2]
            Wn_norm = [w / nyquist for w in Wn]
            return signal.butter(self.order, Wn_norm, btype=self.filter_type,
                                output='sos')

        return None


def sos_group_delay(sos, w, fs=1.0):
    """Calculate group delay for SOS filter"""
    try:
        # 计算频率响应
        _, h = signal.sosfreqz(sos, worN=w, fs=fs)
        phase = np.unwrap(np.angle(h))

        if len(w) > 1:
            # 计算频率差
            dw = np.diff(w)

            # 计算群延迟（使用中心差分或前向差分）
            # 方法：使用中点频率的群延迟
            gd_mid = -np.diff(phase) / dw

            # 插值到原始频率点
            # w_mid 是中点频率
            w_mid = (w[:-1] + w[1:]) / 2

            # 插值回原始频率点
            gd = np.interp(w, w_mid, gd_mid)

            # 处理边界点（使用最近邻）
            if len(gd_mid) > 0:
                gd[0] = gd_mid[0]  # 第一个点使用第一个中点值
                gd[-1] = gd_mid[-1]  # 最后一个点使用最后一个中点值
        else:
            gd = np.zeros_like(w)

        # 转换为样本数（弧度/秒 -> 秒 -> 样本数）
        # gd 当前单位是 弧度/弧度/秒？实际是 -dφ/dω，单位是 秒
        # 需要转换为样本数：乘以 fs 得到样本数
        gd = gd / (2 * np.pi) * fs

        return gd

    except Exception as e:
        print(f"  Warning: Could not compute group delay: {e}")
        return np.zeros_like(w)


def compare_filters(reader, scipy_filter):
    dsp_sos = reader.get_sos()
    scipy_sos = scipy_filter.sos

    print("="*80)
    print("FILTER COMPARISON ANALYSIS")
    print("="*80)

    # 1. Basic Information
    print("\n📊 BASIC INFORMATION")
    print("-"*80)
    print(f"Filter Type:     {reader.filter_type}")
    print(f"Order:           {reader.filter_info.get('order', 'N/A')}")
    print(f"Sample Rate:     {reader.filter_info.get('sample_rate', 'N/A')} Hz")

    if 'cutoff' in reader.filter_info:
        print(f"Cutoff:          {reader.filter_info['cutoff']} Hz")
    if 'center' in reader.filter_info:
        print(f"Center:          {reader.filter_info['center']} Hz")
    if 'bandwidth' in reader.filter_info:
        print(f"Bandwidth:       {reader.filter_info['bandwidth']} Hz")

    print(f"\nSOS Stages:      {len(dsp_sos)} (DspFilters) vs {len(scipy_sos)} (SciPy)")

    # 2. Coefficient Analysis
    print("\n📐 COEFFICIENT ANALYSIS")
    print("-"*80)

    if len(dsp_sos) == len(scipy_sos):
        max_b_diff = 0
        max_a_diff = 0
        mean_b_diff = 0
        mean_a_diff = 0

        for i in range(len(dsp_sos)):
            b_diff = np.abs(dsp_sos[i][:3] - scipy_sos[i][:3])
            a_diff = np.abs(dsp_sos[i][3:] - scipy_sos[i][3:])

            max_b_diff = max(max_b_diff, np.max(b_diff))
            max_a_diff = max(max_a_diff, np.max(a_diff))
            mean_b_diff += np.mean(b_diff)
            mean_a_diff += np.mean(a_diff)

            print(f"\n  Stage {i+1}:")
            print(f"    SciPy:      b=[{scipy_sos[i][0]:>12.8f}, {scipy_sos[i][1]:>12.8f}, {scipy_sos[i][2]:>12.8f}]")
            print(f"                a=[{scipy_sos[i][3]:>12.8f}, {scipy_sos[i][4]:>12.8f}, {scipy_sos[i][5]:>12.8f}]")
            print(f"    DspFilters: b=[{dsp_sos[i][0]:>12.8f}, {dsp_sos[i][1]:>12.8f}, {dsp_sos[i][2]:>12.8f}]")
            print(f"                a=[{dsp_sos[i][3]:>12.8f}, {dsp_sos[i][4]:>12.8f}, {dsp_sos[i][5]:>12.8f}]")
            print(f"    Max Diff:   b={np.max(b_diff):.2e}, a={np.max(a_diff):.2e}")

        mean_b_diff /= len(dsp_sos)
        mean_a_diff /= len(dsp_sos)

        print(f"\n  Summary:")
        print(f"    Mean Coefficient Diff: b={mean_b_diff:.2e}, a={mean_a_diff:.2e}")
        print(f"    Max Coefficient Diff:  b={max_b_diff:.2e}, a={max_a_diff:.2e}")

    # 3. Pole-Zero Analysis
    print("\n🎯 POLE-ZERO ANALYSIS")
    print("-"*80)

    try:
        z_dsp, p_dsp, k_dsp = sos2zpk(dsp_sos)
        z_scipy, p_scipy, k_scipy = sos2zpk(scipy_sos)

        # Sort poles and zeros by magnitude for comparison
        z_dsp_sorted = sorted(z_dsp, key=lambda x: np.abs(x))
        p_dsp_sorted = sorted(p_dsp, key=lambda x: np.abs(x))
        z_scipy_sorted = sorted(z_scipy, key=lambda x: np.abs(x))
        p_scipy_sorted = sorted(p_scipy, key=lambda x: np.abs(x))

        print(f"\n  Poles:")
        for i, (p_d, p_s) in enumerate(zip(p_dsp_sorted[:5], p_scipy_sorted[:5])):
            if i < 5:  # Show first 5 poles
                print(f"    Pole {i+1}: DspFilters={p_d.real:>8.4f}{p_d.imag:>+8.4f}j, "
                      f"SciPy={p_s.real:>8.4f}{p_s.imag:>+8.4f}j")

        # Pole magnitude and angle analysis
        p_mag_dsp = np.abs(p_dsp)
        p_mag_scipy = np.abs(p_scipy)
        p_angle_dsp = np.angle(p_dsp)
        p_angle_scipy = np.angle(p_scipy)

        print(f"\n  Pole Magnitudes:")
        print(f"    Max:    {np.max(p_mag_dsp):.6f} (DspFilters) vs {np.max(p_mag_scipy):.6f} (SciPy)")
        print(f"    Min:    {np.min(p_mag_dsp):.6f} vs {np.min(p_mag_scipy):.6f}")
        print(f"    Mean:   {np.mean(p_mag_dsp):.6f} vs {np.mean(p_mag_scipy):.6f}")
        print(f"    StdDev: {np.std(p_mag_dsp):.6f} vs {np.std(p_mag_scipy):.6f}")

        # Stability check
        print(f"\n  Stability Check:")
        dsp_stable = np.all(np.abs(p_dsp) < 1)
        scipy_stable = np.all(np.abs(p_scipy) < 1)
        print(f"    DspFilters: {'✓ STABLE' if dsp_stable else '✗ UNSTABLE'}")
        print(f"    SciPy:      {'✓ STABLE' if scipy_stable else '✗ UNSTABLE'}")

    except Exception as e:
        print(f"  Could not perform pole-zero analysis: {e}")

    # 4. Frequency Response Analysis
    print("\n📈 FREQUENCY RESPONSE ANALYSIS")
    print("-"*80)

    # Use higher resolution for accurate comparison
    n_points = 4096
    fs = reader.filter_info.get('sample_rate', 44100)

    w_dsp, h_dsp = signal.sosfreqz(dsp_sos, worN=n_points, fs=fs)
    w_scipy, h_scipy = signal.sosfreqz(scipy_sos, worN=n_points, fs=fs)

    mag_dsp = 20 * np.log10(np.abs(h_dsp) + 1e-12)
    mag_scipy = 20 * np.log10(np.abs(h_scipy) + 1e-12)
    phase_dsp = np.unwrap(np.angle(h_dsp)) * 180 / np.pi
    phase_scipy = np.unwrap(np.angle(h_scipy)) * 180 / np.pi

    # Magnitude error
    mag_error = np.abs(mag_dsp - mag_scipy)
    mean_mag_error = np.mean(mag_error)
    max_mag_error = np.max(mag_error)
    rms_mag_error = np.sqrt(np.mean(mag_error**2))

    print(f"\n  Magnitude Response Error:")
    print(f"    Mean Error:  {mean_mag_error:.6f} dB")
    print(f"    RMS Error:   {rms_mag_error:.6f} dB")
    print(f"    Max Error:   {max_mag_error:.6f} dB")

    # Phase error
    phase_error = np.abs(phase_dsp - phase_scipy)
    mean_phase_error = np.mean(phase_error)
    max_phase_error = np.max(phase_error)

    print(f"\n  Phase Response Error:")
    print(f"    Mean Error:  {mean_phase_error:.4f} degrees")
    print(f"    Max Error:   {max_phase_error:.4f} degrees")

    # Group delay
    gd_dsp = sos_group_delay(dsp_sos, w_dsp, fs)
    gd_scipy = sos_group_delay(scipy_sos, w_scipy, fs)
    gd_error = np.abs(gd_dsp - gd_scipy)
    mean_gd_error = np.mean(gd_error)
    max_gd_error = np.max(gd_error)

    print(f"\n  Group Delay Error:")
    print(f"    Mean Error:  {mean_gd_error:.4f} samples")
    print(f"    Max Error:   {max_gd_error:.4f} samples")

    # Frequency band analysis
    if 'cutoff' in reader.filter_info:
        cutoff = reader.filter_info['cutoff']
        # Find indices near cutoff
        idx_cutoff = np.argmin(np.abs(w_dsp - cutoff))
        print(f"\n  At Cutoff Frequency ({cutoff} Hz):")
        print(f"    Magnitude: {mag_dsp[idx_cutoff]:.2f} dB (DspFilters) vs {mag_scipy[idx_cutoff]:.2f} dB (SciPy)")
        print(f"    Phase:     {phase_dsp[idx_cutoff]:.2f}° vs {phase_scipy[idx_cutoff]:.2f}°")

    # 5. Numerical Properties
    print("\n🔢 NUMERICAL PROPERTIES")
    print("-"*80)

    # Check for numerical issues
    dsp_max_coeff = np.max(np.abs(dsp_sos))
    scipy_max_coeff = np.max(np.abs(scipy_sos))

    print(f"\n  Coefficient Range:")
    print(f"    DspFilters: {np.min(dsp_sos):.2e} to {np.max(dsp_sos):.2e}")
    print(f"    SciPy:      {np.min(scipy_sos):.2e} to {np.max(scipy_sos):.2e}")

    # Check for near-unstable poles
    p_mag_dsp = np.abs(p_dsp) if 'p_dsp' in locals() else []
    if len(p_mag_dsp) > 0:
        near_unstable = np.sum(p_mag_dsp > 0.9999)
        if near_unstable > 0:
            print(f"\n  ⚠️  Warning: {near_unstable} poles near unit circle (|p| > 0.9999)")

    # 6. Conclusion
    print("\n📋 CONCLUSION")
    print("-"*80)

    if max_mag_error < 1e-6 and max_phase_error < 1e-4:
        print("✅ Filters are IDENTICAL (within numerical precision)")
        print("   The implementations match perfectly.")
    elif max_mag_error < 1e-3:
        print("✅ Filters are VERY CLOSE")
        print("   Differences are due to numerical rounding and implementation details.")
    elif max_mag_error < 1e-1:
        print("⚠️  Filters show ACCEPTABLE differences")
        print("   May be due to different SOS ordering or numerical methods.")
    else:
        print("❌ Filters show SIGNIFICANT DIFFERENCES")
        print("   Check filter design parameters and SOS ordering.")

    plot_comparison(w_dsp, mag_dsp, phase_dsp, gd_dsp,
            w_scipy, mag_scipy, phase_scipy, gd_scipy,
            mag_error, phase_error)

    print_summary_report(w_dsp, mag_error, phase_error)


def plot_comparison(w_dsp, mag_dsp, phase_dsp, gd_dsp,
                   w_scipy, mag_scipy, phase_scipy, gd_scipy,
                   mag_error, phase_error):
    """
    Comprehensive visualization of filter comparison
    """
    # Create figure with GridSpec for better layout
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)

    # 1. Magnitude Response
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.semilogx(w_dsp, mag_dsp, 'b-', label='DspFilters', linewidth=2, alpha=0.8)
    ax1.semilogx(w_scipy, mag_scipy, 'r--', label='SciPy', linewidth=2, alpha=0.8)
    ax1.set_xlabel('Frequency (Hz)')
    ax1.set_ylabel('Magnitude (dB)')
    ax1.set_title('Magnitude Response Comparison', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='best')
    ax1.set_xlim([min(w_dsp), max(w_dsp)])

    # 2. Phase Response
    ax2 = fig.add_subplot(gs[1, :2])
    ax2.semilogx(w_dsp, phase_dsp, 'b-', label='DspFilters', linewidth=2, alpha=0.8)
    ax2.semilogx(w_scipy, phase_scipy, 'r--', label='SciPy', linewidth=2, alpha=0.8)
    ax2.set_xlabel('Frequency (Hz)')
    ax2.set_ylabel('Phase (degrees)')
    ax2.set_title('Phase Response Comparison', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best')
    ax2.set_xlim([min(w_dsp), max(w_dsp)])

    # 3. Group Delay
    ax3 = fig.add_subplot(gs[2, :2])
    ax3.semilogx(w_dsp, gd_dsp, 'b-', label='DspFilters', linewidth=2, alpha=0.8)
    ax3.semilogx(w_scipy, gd_scipy, 'r--', label='SciPy', linewidth=2, alpha=0.8)
    ax3.set_xlabel('Frequency (Hz)')
    ax3.set_ylabel('Group Delay (samples)')
    ax3.set_title('Group Delay Comparison', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc='best')
    ax3.set_xlim([min(w_dsp), max(w_dsp)])

    # 4. Magnitude Error
    ax4 = fig.add_subplot(gs[0, 2])
    ax4.semilogx(w_dsp, mag_error, 'g-', linewidth=2)
    ax4.set_xlabel('Frequency (Hz)')
    ax4.set_ylabel('Error (dB)')
    ax4.set_title('Magnitude Error', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.set_yscale('log')
    ax4.set_xlim([min(w_dsp), max(w_dsp)])

    # 5. Phase Error
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.semilogx(w_dsp, phase_error, 'orange', linewidth=2)
    ax5.set_xlabel('Frequency (Hz)')
    ax5.set_ylabel('Error (degrees)')
    ax5.set_title('Phase Error', fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3)
    ax5.set_yscale('log')
    ax5.set_xlim([min(w_dsp), max(w_dsp)])

    # 6. Error Statistics
    ax6 = fig.add_subplot(gs[2, 2])
    ax6.axis('off')

    # Add statistics text box
    stats_text = (
        "📊 ERROR STATISTICS\n" +
        "═" * 30 + "\n\n" +
        f"MAGNITUDE ERROR:\n"
        f"  Mean: {np.mean(mag_error):.4f} dB\n"
        f"  RMS:  {np.sqrt(np.mean(mag_error**2)):.4f} dB\n"
        f"  Max:  {np.max(mag_error):.4f} dB\n\n"
        f"PHASE ERROR:\n"
        f"  Mean: {np.mean(phase_error):.4f}°\n"
        f"  RMS:  {np.sqrt(np.mean(phase_error**2)):.4f}°\n"
        f"  Max:  {np.max(phase_error):.4f}°"
    )

    ax6.text(0.05, 0.95, stats_text, transform=ax6.transAxes,
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Add overall title
    fig.suptitle('Digital Filter Implementation Comparison',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout()
    plt.show()

    # Optional: Create additional specialized plots
    create_additional_plots(w_dsp, mag_dsp, mag_scipy, phase_dsp,
                           phase_scipy, mag_error, phase_error)


def create_additional_plots(w_dsp, mag_dsp, mag_scipy, phase_dsp,
                           phase_scipy, mag_error, phase_error):
    """
    Create additional specialized comparison plots
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Zoomed view of magnitude response around transition band
    ax = axes[0, 0]
    # Find the transition region (where magnitude is between -60dB and -3dB)
    mask = (mag_dsp > -60) & (mag_dsp < -3)
    if np.any(mask):
        ax.semilogx(w_dsp[mask], mag_dsp[mask], 'b-', label='DspFilters', linewidth=2)
        ax.semilogx(w_dsp[mask], mag_scipy[mask], 'r--', label='SciPy', linewidth=2)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Magnitude (dB)')
        ax.set_title('Transition Band Detail')
        ax.grid(True, alpha=0.3)
        ax.legend()

    # Passband ripple analysis
    ax = axes[0, 1]
    # Find passband (e.g., below cutoff if lowpass)
    passband_mask = mag_dsp > -3
    if np.any(passband_mask):
        passband_error = mag_error[passband_mask]
        ax.semilogx(w_dsp[passband_mask], passband_error, 'purple', linewidth=2)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Error (dB)')
        ax.set_title('Passband Error Detail')
        ax.grid(True, alpha=0.3)
        ax.set_yscale('log')

    # Phase difference vs frequency
    ax = axes[1, 0]
    ax.semilogx(w_dsp, phase_error, 'orange', linewidth=2)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Phase Difference (degrees)')
    ax.set_title('Phase Difference vs Frequency')
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')

    # Cumulative error
    ax = axes[1, 1]
    cumulative_mag_error = np.cumsum(mag_error) / len(mag_error)
    cumulative_phase_error = np.cumsum(phase_error) / len(phase_error)
    ax.semilogx(w_dsp, cumulative_mag_error, 'green', label='Magnitude', linewidth=2)
    ax.semilogx(w_dsp, cumulative_phase_error, 'red', label='Phase', linewidth=2)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Cumulative Mean Error')
    ax.set_title('Cumulative Error Analysis')
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.show()


def print_summary_report(w_dsp, mag_error, phase_error):
    """
    Print a concise summary report
    """
    print("\n" + "="*80)
    print("SUMMARY REPORT")
    print("="*80)

    # Frequency bands analysis
    bands = [
        ('Full Band', w_dsp[0], w_dsp[-1]),
        ('Low Frequency', w_dsp[0], w_dsp[-1]/10),
        ('Mid Frequency', w_dsp[-1]/10, w_dsp[-1]/2),
        ('High Frequency', w_dsp[-1]/2, w_dsp[-1])
    ]

    for band_name, f_low, f_high in bands:
        mask = (w_dsp >= f_low) & (w_dsp <= f_high)
        if np.any(mask):
            mag_err_band = mag_error[mask]
            phase_err_band = phase_error[mask]
            print(f"\n{band_name}:")
            print(f"  Magnitude Error: mean={np.mean(mag_err_band):.4f} dB, "
                  f"max={np.max(mag_err_band):.4f} dB")
            print(f"  Phase Error:     mean={np.mean(phase_err_band):.4f}°, "
                  f"max={np.max(phase_err_band):.4f}°")


def main():
    parser = argparse.ArgumentParser(description='Compare DspFilters with SciPy Butterworth')
    parser.add_argument('file', help='Filter parameter file')

    args = parser.parse_args()

    if not Path(args.file).exists():
        print(f"Error: File {args.file} does not exist")
        return 1

    reader = DspFiltersReader(args.file)

    if not reader.filter_type:
        print("Error: Cannot determine filter type from file")
        return 1

    params = dict()
    if reader.filter_type in ['lowpass', 'highpass']:
        params['cutoff'] = reader.filter_info.get('cutoff', 1000)
    else:
        params['center'] = reader.filter_info.get('center', 1000)
        params['bandwidth'] = reader.filter_info.get('bandwidth', 100)

    scipy_filter = SciPyButterworth(
        filter_type=reader.filter_type,
        order=reader.filter_info.get('order', 4),
        sample_rate=reader.filter_info.get('sample_rate', 44100),
        **{k: v for k, v in params.items()}
    )

    compare_filters(reader, scipy_filter)

    return 0


if __name__ == "__main__":
    sys.exit(main())
