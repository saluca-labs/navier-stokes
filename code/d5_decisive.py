r"""d5_decisive.py — ONE well-resolved 512^3 Re=3200 run for a trustworthy extreme-tail point.

At 512^3, Re=3200 gives k_max*eta ~ 1.5 (vs the under-resolved 0.86 at 384^3). Compares the extreme-tail
depletion (top 0.01% |omega|) at GOOD resolution against the under-resolved 384^3 value (0.099):
  tail_512 > 0.099  => under-resolution was clipping extremes; true tail DE-depletes (toward singularity).
  tail_512 ~ 0.099  => the fall-back was real (extremes re-deplete; toward regularity).
If NS_AFTER_CMD is set, runs that command when done (the original run used it to restart a
local service that shared the GPU). Run in background.
"""
import os, time, subprocess
import d4_extreme as d4

print("D5 — decisive 512^3 Re=3200 (well-resolved extreme tail)...", flush=True)
t0 = time.time()
try:
    r = d4.run_case(512, 3200, 10.0, 0.0035)
    out = (f"Re={r['Re']:.0f} N={r['N']} k_max*eta={r['kme']:.2f} t_peak={r['tpk']:.2f} enstr_pk={r['enspk']:.3f}\n"
           f"  mean={r['mean']:.3f}  top1%={r['ex'][0]:.3f}  top0.1%={r['ex'][1]:.3f}  top0.01%={r['ex'][2]:.3f}\n"
           f"  vs under-resolved 384^3 Re=3200 tail was 0.099 (k_max*eta 0.86)\n"
           f"  elapsed {time.time()-t0:.0f}s")
    print(out, flush=True)
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "experiments")
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "d5_result.txt"), "w", encoding="utf-8") as f:
        f.write(out + "\n")
except Exception as e:
    print("D5 FAILED:", str(e)[:200], flush=True)
finally:
    after = os.environ.get("NS_AFTER_CMD")
    if after:
        subprocess.Popen(after, shell=True)
        print("NS_AFTER_CMD started.", flush=True)
