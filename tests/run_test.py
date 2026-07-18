import argparse
import subprocess
from pathlib import Path
import sys
import os

if __name__ != "__main__":
    sys.exit()

parser = argparse.ArgumentParser()
parser.add_argument("--exe", required=True)
parser.add_argument("--config", required=True)
parser.add_argument("--plot-script", required=True)

args = parser.parse_args()

exe = Path(args.exe)
config = Path(args.config)
plot_script = Path(args.plot_script)

print(f"Running computation for ContinuumSystem:")
print(f"  executable : {exe}")
print(f"  config     : {config}")

simulation_cmd = [
    str(exe),
    str(config)
]
result = subprocess.run(simulation_cmd)

if result.returncode != 0:
    print("Simulation failed.", file=sys.stderr)
    sys.exit(result.returncode)

print("Simulation finished successfully.")
print("Running plot script.")

plot_cmd = [
    sys.executable,
    str(plot_script),
    os.path.splitext(os.path.basename(str(config)))[0]
]
result = subprocess.run(plot_cmd)

if result.returncode != 0:
    print("Plot script failed.", file=sys.stderr)
    sys.exit(result.returncode)

print("Test workflow finished.")
sys.exit(0)