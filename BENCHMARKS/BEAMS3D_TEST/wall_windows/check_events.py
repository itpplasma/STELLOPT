#!/usr/bin/env python3
"""Independent analytic event and rank/reload oracle for synthetic planes."""
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--faces', type=int, choices=[1, 16], required=True)
    parser.add_argument('--rays', type=Path, required=True)
    parser.add_argument('--ranks', type=int, choices=[1, 2], required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--log', type=Path, required=True)
    args = parser.parse_args()
    rays = args.rays.read_text().splitlines()
    assert int(rays[0]) == args.faces and len(rays) == args.faces + 1
    targets = []
    for ray in rays[1:]:
        values = list(map(float, ray.split()))
        assert len(values) == 6
        targets.append([(a + b) / 2 for a, b in zip(values[:3], values[3:])])
    reference = None
    log = args.log.read_text().splitlines()
    for rank in range(args.ranks):
        for cycle in range(1, 4):
            suffix = f'.rank{rank}' if args.ranks > 1 else ''
            path = args.output.with_name(args.output.name + suffix + f'.cycle{cycle}')
            raw = path.read_bytes()
            rows = raw.decode().splitlines()
            assert list(map(int, rows[0].split())) == [args.faces, 3 * args.faces, args.faces]
            assert len(rows) == args.faces + 1
            for row, (raw_event, target) in enumerate(zip(rows[1:], targets), start=1):
                event = list(map(float, raw_event.split()))
                assert len(event) == 7 and event[:3] == [row, 1, row]
                assert max(abs(x - y) for x, y in zip(event[3:], [*target, .5])) <= 1e-12
            if reference is None:
                reference = raw
            assert raw == reference
            assert log.count(f'WALL_CYCLE {cycle} RANK {rank}') == 1
    print('Analytic event, original face ID, rank identity and three reloads PASS')


if __name__ == '__main__':
    main()
