# Shared candidate-window regression

Accelerated wall blocks retain their original ordered face-ID slices in shared
chunks, rather than owning one persistent MPI window each. This preserves the
text format, candidate membership and native triangle intersection. It also
covers candidate storage made by the native generated-grid path. Blocks release
their views before `wall_free` collectively frees the pool. Non-MPI allocations
and the default owned allocation of public `INIT_BLOCK` calls are unchanged.

The synthetic input has constructed planes and deliberately broad candidate
lists. It contains no machine geometry, field, orbit or collaborator data.
The two-bin control checks one plane. The 100,000-bin cases use one or sixteen
planes, with the latter storing 1,600,000 candidate IDs across a chunk boundary.
Each plane has an isolated analytic ray, including planes in bins 65,535 and
65,536 on opposite sides of that boundary. Expected original face IDs are 1–16,
physical hits are the ray midpoints, and segment fractions are 0.5.

Build `native_wall_probe.f90` against the current native MPI LIBSTELL using the
project's configured compiler and native link libraries. Run each case at one
and two MPI ranks inside a matching allocation, with three load/collide/free
cycles. Preserve its separate rank/cycle event files. For example, after
placing the compiled driver in an absolute path:

```sh
python3 prepare_plane.py --bins 100000 --faces 16 --output chunks
cd chunks
mpirun -np 2 /absolute/path/wall_probe.x wall.dat native.events 3 > native.log
python3 /absolute/path/check_events.py --faces 16 --ranks 2 \
    --rays rays.dat --output native.events --log native.log
```

Capture the driver's stdout in `native.log`. `check_events.py` independently
checks every face ID, midpoint and fraction within 1e-12, byte identity across
ranks/reloads, and completion of each cleanup cycle. The native non-MPI wall
module compile and existing `check_event_modes.py --case all` are also required.
The latter includes the published external-wall case and its original field
validity and event oracles, independently of this small reader probe.

The diagnostic profiler is specific to matching GNU/Open MPI headers and ABI:

```sh
mpicc -O2 -fPIC -shared trace_windows.c -ldl -o trace_windows.so
```

Apply `LD_PRELOAD` to the native driver, or export it to MPI application ranks
with `mpirun -x LD_PRELOAD=/absolute/path/trace_windows.so`. It transparently
forwards the actual seven-argument `mpi_win_allocate_shared_cptr_` binding and
communicator duplication, recording allocation/CID state and returned failures.
It changes no MPI provider, limit, error policy, allocator or collision backend.
Record the actual linked/loaded providers and matching internal-header identities,
and use a memory guard and forced timeout appropriate to the allocation.

With Open MPI 4.1.4, the unchanged reader's 100,000-bin synthetic diagnostic
fails at allocation attempt 131,050 with communicator position 65,535 equal to
the PML context-ID bound, and communicator duplication returns 17. Both analytic
controls pass with identical output. This demonstrates context-ID exhaustion
for this public reproducer, without identifying an unrecorded failure in any
other application. The chunked candidate reader passes the one/two-rank,
three-reload matrix and all sixteen isolated ID/boundary-view checks. The
existing event oracles also pass across the retained first three cases and an
identical last-case continuation after an undersized suite timeout. That timeout
is retained rather than reported as a successful full-suite invocation.

The repository has no `INIT_BLOCK` or `BLOCK_DESTROY` callers outside `wall_mod`.
Their public source interface remains compatible: only the internal reader opts
into pooled views, and external calls default to their original owned allocation.
Original MPI allocation helpers and fatal error behavior remain unchanged.
