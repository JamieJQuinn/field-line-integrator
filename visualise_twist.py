import matplotlib.pyplot as plt
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(description='Visualise twist data')
    parser.add_argument('inputs', nargs='*', help='Input npy files')
    args = parser.parse_args()

    twist_grids = [np.load(filename) for filename in args.inputs]

    n_columns = 1
    n_rows = int(len(twist_grids) / n_columns) + (len(twist_grids)%n_columns)

    # fig, axes = plt.subplots(n_rows, n_columns, sharex=True, sharey=True, figsize=(2*n_columns, 2*n_rows))
    # axes = axes.flatten()
    for filename,twist in zip(args.inputs, twist_grids):
        plt.imshow(twist.transpose())
        plt.tight_layout()
        plt.savefig(filename + ".png")

        # plt.show()

main()
