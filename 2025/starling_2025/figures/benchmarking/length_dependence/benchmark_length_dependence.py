import argparse
import time
from starling import generate
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Run the Starling generate function with adjustable parameters.")
    parser.add_argument("--device", type=str, default="mps", help="Device to use for computation (e.g., cpu, cuda, mps)")
    parser.add_argument("--sleepytime", type=int, default=120, help="Time to sleep between runs in seconds")
    args = parser.parse_args()

    block = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI"
    all_times = []
    lengths = []

    # warm up -  run this to make sure everything is loaded
    generate(block, device=args.device)
    
    for r in range(1, 8):
        print(f"On {r} of 7")
        s = block * r
        lengths.append(len(s))

        start = time.time()
        
        # Generate ensemble using default settings
        generate(s, device=args.device)
        
        end = time.time()
        all_times.append(end - start)
                
        time.sleep(args.sleepytime)  # Sleep to avoid rate limiting

    np.savetxt(f'time_vs_seqlen_{args.device}.csv', np.array([lengths, all_times]).T)
    print("Done")

if __name__ == "__main__":
    main()