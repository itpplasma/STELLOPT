#!/usr/bin/env python3
"""Write a synthetic one-triangle accelerated wall with deliberately broad bins."""
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bins', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--faces', type=int, choices=[1, 16], default=1)
    args = parser.parse_args()
    assert 1 <= args.bins <= 1_000_000
    if args.faces == 1:
        bins = [0]
    else:
        assert args.bins > 65536
        bins = [round(i * (args.bins - 1) / 13) for i in range(14)]
        bins += [65535, 65536]
        assert len(set(bins)) == 16
    # The final two analytic hits exercise slices on both sides of a chunk.
    positions = [index + .5 for index in bins]
    args.output.mkdir()
    with (args.output / 'wall.dat').open('x') as stream:
        stream.write('Synthetic planes\nNative accelerated text\n0 0\n')
        stream.write(f'{3 * args.faces} {args.faces}\n')
        for face in range(args.faces):
            x = positions[face]
            stream.write(f'{x:.2f} 0 0\n{x:.2f} 1 0\n{x:.2f} 0 1\n')
        for face in range(args.faces):
            stream.write(f'{3 * face + 1} {3 * face + 2} {3 * face + 3}\n')
        stream.write(f'{args.bins} 1 1 1\n1 {args.bins} 1 1\n')
        for index in range(args.bins):
            stream.write(f'{index} {index + 1} 0 1 0 1\n{args.faces}\n')
            for face in range(1, args.faces + 1):
                stream.write(f'{face}\n')
    with (args.output / 'rays.dat').open('x') as stream:
        stream.write(f'{args.faces}\n')
        for x in positions:
            stream.write(f'{x - .25} 0.2 0.2 {x + .25} 0.2 0.2\n')


if __name__ == '__main__':
    main()
