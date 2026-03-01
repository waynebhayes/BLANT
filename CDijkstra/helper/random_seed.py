import os
import random

# Directories
SIM_DIR = "/extra/wayne1/src/bionets/SANA.github/sims/BioGRID-3.0.064/importance"
SEED_DIR = "/home/yiranj12/BLANT-workspace/BLANT/CDijkstra/test/importancesim_seed/100_random_seed"
os.makedirs(SEED_DIR, exist_ok=True)
Delta = 0.5

# Iterate over all .sim files
for filename in os.listdir(SIM_DIR):
    if not filename.endswith(".sim"):
        continue

    sim_path = os.path.join(SIM_DIR, filename)
    base_name = filename.replace(".sim", "")

    # get max from the top (first valid) line: third token
    max = None
    with open(sim_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    max = float(parts[-1])
                    break
                except ValueError:
                    continue
    if max is None:
        print(f"{filename}: no valid lines; skipping")
        continue

    upper_bound = max
    bottom_bound = max - Delta

    for i in range(1, 101):
        threshold = round(random.uniform(bottom_bound, upper_bound), 6)

        subdir = os.path.join(SEED_DIR, base_name)
        os.makedirs(subdir, exist_ok=True)
        out_path = os.path.join(subdir, f"{base_name}-{i}.txt")

        prev_pair = None
        selected = None

        with open(sim_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue

                node1, node2, score_str = parts[0], parts[1], parts[-1]
                try:
                    score = float(score_str)
                except ValueError:
                    continue

                if score == threshold:
                    selected = (node1, node2, score)
                    break
                elif score < threshold:
                    selected = prev_pair  # pick previous line
                    break

                prev_pair = (node1, node2, score)

        if selected:
            with open(out_path, 'w') as out:
                out.write(f"{selected[0]}\t{selected[1]}\n")
            if selected[2] == threshold:
                print(f"{filename} (times{i}): exact match {selected[0]}-{selected[1]} score={selected[2]}")
            else:
                print(f"{filename} (times{i}): fallback {selected[0]}-{selected[1]} score={selected[2]} (threshold={threshold})")
        else:
            print(f"{filename} (times{i}): no suitable seed found for threshold {threshold}")

