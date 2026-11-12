# Paper 8 patch to the ntn-v2x module (Nov 2026)

`ns-3-dev` is shared with running campaigns, so this revision did not touch it or its build. All of the paper's ns-3 data come
from one standalone binary, `../bin/p8-canyon-relay`, built by `build.sh` against the optimized libraries in
`ns-3-dev/build/lib`. Its SHA-256 prefix is recorded per row in `../out/ns3_canyon5.jsonl`.

- `orig/`: the example and sidelink model as copied from `contrib/ntn-v2x` on 4 Nov 2026.
- `apply_edits.py`: produces the patched copies from `orig/`; `p8.patch` is the resulting diff.
- New flags: `--ulKbps --grantRtt --ulMaxQueueS --relaySetupMs --nbAgeS --kappa --mobTrace --facadeFile --gapEveryM --gapLenM
  --slUrban --slNlosP --alignResel`; new policies `rand2 rand3 roundrobin hybrid anyfwd anyp`.

The model classes are compiled under other names (`-DNtnSlUeMac=NtnSlUeMacP8 -DNtnSlChannel=NtnSlChannelP8`, with `P8_STANDALONE`
renaming the TypeIds and log component), so the library keeps its own classes and nothing overrides them.

Regression: with every new flag at its default, the binary reproduced the stored rows of the old binary exactly (urban, rep 0: `random`
at 100 ms, 188737 generated, 127506 delivered, 12801 relayed; `multisat` at 1 s, 19846 generated, 12066 delivered).

**To merge at the shared build:** apply `p8.patch` under `contrib/ntn-v2x` (model and example paths), drop the `P8_STANDALONE`
blocks and change the example's include back to `"ns3/ntn-nr-sidelink.h"`, then rebuild.
