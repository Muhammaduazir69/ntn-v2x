#!/usr/bin/env python3
"""Round 4 edits on the round-3 example (build 4):
1. The neighbors a message names travel in its 32-bit tag (up to three ids of 10 bits) for relayk and relayrand with
   relayK <= 3, as rand2/rand3 already did, instead of being read from the source's state at decode. For more than three
   names (randk with 5 or 8) the list is still read from the source's state.
2. Two announcer rules the model suggests: relayfar (the farthest announcing neighbor in range) and relaylb (the nearest
   announcing neighbor at least one building length away along the street, else the nearest announcer)."""
import os
H = os.path.dirname(os.path.abspath(__file__))
s = open(f"{H}/orig/ntn-v2x-canyon-relay.cc").read()
def rep(a, b):
    global s
    assert s.count(a) == 1, a[:70]
    s = s.replace(a, b)
rep('if ((g_policy == "relay" || g_policy == "relayk" || g_policy == "relayrand") && !n.up)',
    'if ((g_policy == "relay" || g_policy == "relayk" || g_policy == "relayrand" || g_policy == "relayfar" || g_policy == "relaylb") && !n.up)')
rep('''        if (d < bestD)
        {
            bestD = d;
            pick = j;
        }''', '''        if (d < bestD)
        {
            bestD = d;
            pick = j;
        }
        if (d > worstD)
        {
            worstD = d;
            pickFarthest = j;
        }''')
rep('    double bestD = 1e18, bestFar = 1e18, bestLb = 1e18, bestQ = 2.0, bestQD = 1e18;',
    '    double bestD = 1e18, bestFar = 1e18, bestLb = 1e18, bestQ = 2.0, bestQD = 1e18, worstD = -1.0;\n    uint32_t pickFarthest = kNoRelay;')
rep('''    if (g_policy == "relayk" && !cand.empty())
    {''', '''    if (g_policy == "relayfar" && pickFarthest != kNoRelay)
    {
        return pickFarthest + 1;
    }
    if (g_policy == "relaylb")
    {
        return pickLb != kNoRelay ? pickLb + 1 : (pick == kNoRelay ? kNoRelay : pick + 1);
    }
    if (g_policy == "relayk" && !cand.empty())
    {''')
# names in the tag: the same draws as before (same generator, same order), only the carrier changes
rep('''            nm.push_back(byD[m].second);
        }
        return kList;''', '''            nm.push_back(byD[m].second);
        }
        if (g_relayK <= 3)
        {
            uint32_t tagv = 0;
            for (uint32_t m = 0; m < nm.size(); ++m)
            {
                tagv |= (nm[m] + 1) << (10 * m);
            }
            nm.clear();
            return tagv;
        }
        return kList;''')
rep('''            nm.push_back(cand[idx]);
            cand.erase(cand.begin() + idx);
        }
        return kList;''', '''            nm.push_back(cand[idx]);
            cand.erase(cand.begin() + idx);
        }
        if (g_relayK <= 3)
        {
            uint32_t tagv = 0;
            for (uint32_t m = 0; m < nm.size(); ++m)
            {
                tagv |= (nm[m] + 1) << (10 * m);
            }
            nm.clear();
            return tagv;
        }
        return kList;''')
rep('g_policy == "relayrand" || g_policy == "randk";', 'g_policy == "relayrand" || g_policy == "randk" ||\n                          g_policy == "relayfar" || g_policy == "relaylb";')
open(f"{H}/ntn-v2x-canyon-relay.cc", "w").write(s)
