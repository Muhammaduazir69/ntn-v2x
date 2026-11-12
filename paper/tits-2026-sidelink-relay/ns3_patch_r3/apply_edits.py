#!/usr/bin/env python3
"""Round 3 edit: record the slant range on the vehicle, so that --grantRtt takes effect. In the earlier binaries v.slant stayed 0,
which made the scheduling-request round trip zero whatever the flag said."""
import os
H = os.path.dirname(os.path.abspath(__file__))
s = open(f"{H}/orig/ntn-v2x-canyon-relay.cc").read()
a = '            v.chan->SetAttribute("Delay", TimeValue(Seconds(2 * slant / kC)));\n'
assert s.count(a) == 1
s = s.replace(a, a + "            v.slant = slant;\n")
open(f"{H}/ntn-v2x-canyon-relay.cc", "w").write(s)
