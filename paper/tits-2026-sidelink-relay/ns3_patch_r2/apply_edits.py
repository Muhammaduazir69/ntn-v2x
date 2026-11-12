#!/usr/bin/env python3
"""Round 2 edits on top of the first patched example (orig/ = ns3_patch/ntn-v2x-canyon-relay.cc). New flags default to the old
behavior: --relayK, --ulPrio, --vMin, --vMax; new policies relayk, relayrand, randk, anypn."""
import os
H = os.path.dirname(os.path.abspath(__file__))
s = open(f"{H}/orig/ntn-v2x-canyon-relay.cc").read()
def rep(a, b):
    global s
    assert s.count(a) == 1, a[:70]
    s = s.replace(a, b)
rep("constexpr uint32_t kAny = 0xFFFFFFFE;\n", "constexpr uint32_t kAny = 0xFFFFFFFE;\nconstexpr uint32_t kList = 0xFFFFFFFD; // the broadcast names a list of neighbors (g_named)\n")
rep("double g_traceDt = 1.0;\n", """double g_traceDt = 1.0;
uint32_t g_relayK = 1;     // relayk, relayrand, randk: how many neighbors a message names
bool g_ulPrio = false;     // uplink discipline: a vehicle's own messages before forwards, instead of first come first served
double g_vMin = 6.0, g_vMax = 14.0; // constant-speed range (m/s)
std::map<uint32_t, std::vector<uint32_t>> g_named; // per source, the neighbors named by its current message
struct UlItem
{
    uint32_t i;
    Ptr<Packet> p;
    uint16_t port;
    double t0;
};
std::map<int, std::array<std::deque<UlItem>, 2>> g_ulQ; // per satellite: own messages, forwards
std::map<int, bool> g_ulBusy;
""")
rep("#include <algorithm>\n", "#include <algorithm>\n#include <array>\n#include <deque>\n")
# priority discipline
rep("""void
UlSend(uint32_t i, Ptr<Packet> p, uint16_t port)
{""", """void UlServe(int sat);

void
UlDone(uint32_t i, Ptr<Packet> p, uint16_t port, double t0, int sat)
{
    UlFinish(i, p, port, t0);
    UlServe(sat);
}

/// Strict priority, not preemptive: own messages first, each class in arrival order; a message that has waited too long is dropped.
void
UlServe(int sat)
{
    auto& q = g_ulQ[sat];
    const double now = Simulator::Now().GetSeconds();
    for (auto& cls : q)
    {
        while (!cls.empty())
        {
            UlItem it = cls.front();
            cls.pop_front();
            if (now - it.t0 > g_ulMaxQ)
            {
                ++g_ulDropQ;
                continue;
            }
            g_ulBusy[sat] = true;
            const double tx = (it.p->GetSize() + 28) * 8.0 / (g_ulKbps * 1e3);
            Simulator::Schedule(Seconds(tx), &UlDone, it.i, it.p, it.port, it.t0, sat);
            return;
        }
    }
    g_ulBusy[sat] = false;
}

void
UlEnqueue(uint32_t i, Ptr<Packet> p, uint16_t port, double t0)
{
    const int sat = g_v[i].sat;
    g_ulQ[sat][port == kPort ? 0 : 1].push_back({i, p, port, t0});
    if (!g_ulBusy[sat])
    {
        UlServe(sat);
    }
}

void
UlSend(uint32_t i, Ptr<Packet> p, uint16_t port)
{""")
rep("""    double t = now + g_grantRtt * rtt;
    if (g_ulKbps > 0)
    {""", """    double t = now + g_grantRtt * rtt;
    if (g_ulPrio && g_ulKbps > 0)
    {
        Simulator::Schedule(Seconds(t - now), &UlEnqueue, i, p, port, now);
        return;
    }
    if (g_ulKbps > 0)
    {""")
# relay choice
rep("""        if (g_policy == "relay" && !n.up)
        {""", """        if ((g_policy == "relay" || g_policy == "relayk" || g_policy == "relayrand") && !n.up)
        {""")
rep("""        cand.push_back(j);
""", """        cand.push_back(j);
        byD.emplace_back(d, j);
""")
rep("""    std::vector<uint32_t> cand;
""", """    std::vector<uint32_t> cand;
    std::vector<std::pair<double, uint32_t>> byD;
""")
rep("""    if ((g_policy == "rand2" || g_policy == "rand3") && !cand.empty())
    {""", """    if (g_policy == "relayk" && !cand.empty())
    {
        std::sort(byD.begin(), byD.end());
        auto& nm = g_named[src];
        nm.clear();
        for (uint32_t m = 0; m < g_relayK && m < byD.size(); ++m)
        {
            nm.push_back(byD[m].second);
        }
        return kList;
    }
    if ((g_policy == "relayrand" || g_policy == "randk") && !cand.empty())
    {
        static Ptr<UniformRandomVariable> ul = CreateObject<UniformRandomVariable>();
        auto& nm = g_named[src];
        nm.clear();
        for (uint32_t m = 0; m < g_relayK && !cand.empty(); ++m)
        {
            const uint32_t idx = ul->GetInteger(0, cand.size() - 1);
            nm.push_back(cand[idx]);
            cand.erase(cand.begin() + idx);
        }
        return kList;
    }
    if ((g_policy == "rand2" || g_policy == "rand3") && !cand.empty())
    {""")
rep("""                          g_policy == "roundrobin" || g_policy == "hybrid";""", """                          g_policy == "roundrobin" || g_policy == "hybrid" || g_policy == "relayk" ||
                          g_policy == "relayrand" || g_policy == "randk";""")
rep("""    else if (!v.up && (g_policy == "anyfwd" || g_policy == "anyp"))""", """    else if (!v.up && (g_policy == "anyfwd" || g_policy == "anyp" || g_policy == "anypn"))""")
rep("""        if (g_policy == "anyp")
        {""", """        if (g_policy == "anyp" || g_policy == "anypn")
        {""")
rep("""                ups += n.up && now - n.t <= g_nbAge""", """                ups += (n.up || g_policy == "anypn") && now - n.t <= g_nbAge""")
rep("""    else if (tg != kNoRelay)
    {""", """    else if (tg == kList)
    {
        const auto& nm = g_named[from];
        named = std::find(nm.begin(), nm.end(), self) != nm.end();
    }
    else if (tg != kNoRelay)
    {""")
rep("""    cmd.AddValue("kappa", "anyp: expected number of forwards", g_kappa);
""", """    cmd.AddValue("kappa", "anyp: expected number of forwards", g_kappa);
    cmd.AddValue("relayK", "relayk, relayrand, randk: neighbors named per message", g_relayK);
    cmd.AddValue("ulPrio", "own uplink messages before forwards", g_ulPrio);
    cmd.AddValue("vMin", "constant speed range, low (m/s)", g_vMin);
    cmd.AddValue("vMax", "constant speed range, high (m/s)", g_vMax);
""")
rep("""        v.v = u->GetValue(6, 14) * (v.y < g_width / 2 ? 1 : -1);""", """        v.v = u->GetValue(g_vMin, g_vMax) * (v.y < g_width / 2 ? 1 : -1);""")
open(f"{H}/ntn-v2x-canyon-relay.cc", "w").write(s)
