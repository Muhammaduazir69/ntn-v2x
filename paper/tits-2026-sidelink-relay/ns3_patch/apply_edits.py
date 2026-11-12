#!/usr/bin/env python3
"""Produce the patched copies of the sidelink model and the canyon example from orig/. Every edit must match once."""
import os
H = os.path.dirname(os.path.abspath(__file__))
def rep(s, a, b, n=1):
    assert s.count(a) >= 1 and (n == 0 or s.count(a) == n), a[:70]
    return s.replace(a, b)

# ---------- sidelink header
s = open(f"{H}/orig/ntn-nr-sidelink.h").read()
s = rep(s, "    NtnSlGrant SelectResource(uint64_t nowSlot);",
"""    NtnSlGrant SelectResource(uint64_t nowSlot);
    /// P8: selection window after the arrival slot of the next periodic message; sensing is pruned at nowSlot.
    NtnSlGrant SelectResourceAt(uint64_t nowSlot, uint64_t arrivalSlot);""")
s = rep(s, "    void SetRsrpThresholdDbm(double dbm) { m_rsrpThreshDbm = dbm; }",
"""    void SetRsrpThresholdDbm(double dbm) { m_rsrpThreshDbm = dbm; }
    /// P8: message arrivals on a fixed period with a random phase, so reselection does not fire a message early.
    void SetAlignReselection(bool on) { m_align = on; }""")
s = rep(s, "    uint32_t m_pktBytes{190};\n",
"""    uint32_t m_pktBytes{190};
    bool m_align{false};        //!< P8: align reselection with the message period
    uint64_t m_arrival{0};      //!< P8: arrival slot of the current message
""")
s = rep(s, "    void SetDecodeSinrDb(double db) { m_decodeSinrDb = db; }",
"""    void SetDecodeSinrDb(double db) { m_decodeSinrDb = db; }
    /// P8: TR 37.885 urban channel (LOS, NLOSv, NLOS per pair of vehicles) instead of the highway LOS law.
    void SetUrbanChannel(bool on, double nlosP = 0.0);""")
s = rep(s, "    Ptr<NormalRandomVariable> m_nlosvRv;\n",
"""    Ptr<NormalRandomVariable> m_nlosvRv;
    // P8: urban channel; per pair of mobility models, the condition and the distance at which it was drawn. It is
    // kept until the pair's distance has changed by a shadowing correlation distance (TR 37.885: 10 to 13 m).
    bool m_urban{false};
    double m_nlosP{0.0};
    struct PairState
    {
        int cond;
        double d0;
    };
    mutable std::map<std::pair<const void*, const void*>, PairState> m_pair;
    Ptr<UniformRandomVariable> m_condRv;
    Ptr<NormalRandomVariable> m_stdNormal;
    double UrbanLossDb(Ptr<MobilityModel> tx, Ptr<MobilityModel> rx) const;
""")
open(f"{H}/ntn-nr-sidelink.h", "w").write(s)

# ---------- sidelink model
s = open(f"{H}/orig/ntn-nr-sidelink.cc").read()
s = rep(s, 'NS_LOG_COMPONENT_DEFINE("NtnNrSidelink");',
'''#ifdef P8_STANDALONE
NS_LOG_COMPONENT_DEFINE("NtnNrSidelinkP8");
#else
NS_LOG_COMPONENT_DEFINE("NtnNrSidelink");
#endif''')
s = rep(s, 'TypeId("ns3::NtnSlUeMac")', 'TypeId(P8_NAME("ns3::NtnSlUeMac"))')
s = rep(s, 'TypeId("ns3::NtnSlChannel")', 'TypeId(P8_NAME("ns3::NtnSlChannel"))')
s = rep(s, "namespace ns3\n{\n", '''#ifdef P8_STANDALONE
#define P8_NAME(x) x "P8"
#else
#define P8_NAME(x) x
#endif

namespace ns3
{
''')
s = rep(s, """NtnSlGrant
NtnSlUeMac::SelectResource(uint64_t nowSlot)
{
    PruneSensing(nowSlot);
    ++m_reselections;

    const uint32_t C = m_pool.ResourcesPerSlot();
    const uint64_t winStart = nowSlot + m_t1Slots;
    const uint64_t winEnd = nowSlot + m_t2Slots; // inclusive
""", """NtnSlGrant
NtnSlUeMac::SelectResource(uint64_t nowSlot)
{
    return SelectResourceAt(nowSlot, nowSlot);
}

NtnSlGrant
NtnSlUeMac::SelectResourceAt(uint64_t nowSlot, uint64_t arrivalSlot)
{
    PruneSensing(nowSlot);
    ++m_reselections;

    const uint32_t C = m_pool.ResourcesPerSlot();
    const uint64_t winStart = arrivalSlot + m_t1Slots;
    // the packet delay budget of a periodic message does not exceed its period (the aligned mode only)
    const uint64_t winEnd = arrivalSlot + (m_align ? std::min(m_t2Slots, m_prsvpSlots) : m_t2Slots); // inclusive
""")
s = rep(s, """double
NtnSlChannel::SlRsrpDbm(Ptr<MobilityModel> tx, Ptr<MobilityModel> rx) const
{
    if (m_useTr37885)
""", """void
NtnSlChannel::SetUrbanChannel(bool on, double nlosP)
{
    m_urban = on;
    m_nlosP = nlosP;
    if (on && !m_condRv)
    {
        m_condRv = CreateObject<UniformRandomVariable>();
        m_stdNormal = CreateObject<NormalRandomVariable>();
    }
}

double
NtnSlChannel::UrbanLossDb(Ptr<MobilityModel> tx, Ptr<MobilityModel> rx) const
{
    // TR 37.885 urban: Table 6.2-1 gives P(LOS) = min(1, 1.05 exp(-0.0114 d)) for a link no building obstructs, the
    // remainder being NLOSv (a blocking vehicle), as ns-3's ThreeGppV2vUrbanChannelConditionModel implements; two
    // vehicles on one straight street are such a link. Table 6.2.1-1 gives the path loss and shadowing. m_nlosP is a
    // stress knob: the share of pairs given the NLOS law regardless of the geometry.
    const double d = std::max(tx->GetDistanceFrom(rx), 1.0);
    const void* pa = PeekPointer(tx);
    const void* pb = PeekPointer(rx);
    const auto key = pa < pb ? std::make_pair(pa, pb) : std::make_pair(pb, pa);
    auto it = m_pair.find(key);
    if (it == m_pair.end() || std::abs(d - it->second.d0) > 10.0)
    {
        const double pLos = std::min(1.0, 1.05 * std::exp(-0.0114 * d));
        const double stress = m_condRv->GetValue();
        const int cond = stress < m_nlosP ? 2 : (m_condRv->GetValue() < pLos ? 0 : 1);
        m_pair[key] = {cond, d};
        it = m_pair.find(key);
    }
    const double fc = m_carrierHz / 1e9;
    double pl;
    if (it->second.cond == 2)
    {
        pl = 36.85 + 30.0 * std::log10(d) + 18.9 * std::log10(fc) + 4.0 * m_stdNormal->GetValue();
    }
    else
    {
        pl = 38.77 + 16.7 * std::log10(d) + 18.2 * std::log10(fc) + 3.0 * m_stdNormal->GetValue();
        if (it->second.cond == 1)
        {
            // antenna heights equal: case 3 of TR 37.885 clause 6.2.1
            pl += std::max(0.0, 5.0 + std::max(0.0, 15.0 * std::log10(d) - 41.0) + 4.0 * m_stdNormal->GetValue());
        }
    }
    return pl;
}

double
NtnSlChannel::SlRsrpDbm(Ptr<MobilityModel> tx, Ptr<MobilityModel> rx) const
{
    if (m_urban)
    {
        return m_txPowerDbm - UrbanLossDb(tx, rx);
    }
    if (m_useTr37885)
""")
s = rep(s, """        NtnSlGrant g = ue->SelectResource(0);
""", """        NtnSlGrant g;
        if (ue->m_align)
        {
            ue->m_arrival = static_cast<uint64_t>(ue->m_rand->GetInteger(0, ue->m_prsvpSlots - 1));
            g = ue->SelectResourceAt(0, ue->m_arrival);
        }
        else
        {
            g = ue->SelectResource(0);
        }
""")
s = rep(s, """            NtnSlGrant g = t.ue->SelectResource(m_slot);
            m_grants[uid] = g;
""", """            NtnSlGrant g;
            if (t.ue->m_align)
            {
                t.ue->m_arrival += t.ue->m_prsvpSlots;
                g = t.ue->SelectResourceAt(m_slot, t.ue->m_arrival);
            }
            else
            {
                g = t.ue->SelectResource(m_slot);
            }
            m_grants[uid] = g;
""")
# keep the arrival grid moving on kept grants as well
s = rep(s, """            NtnSlGrant g = m_grants[uid];
            g.slot = nextSlot;
            m_grants[uid] = g;
""", """            NtnSlGrant g = m_grants[uid];
            g.slot = nextSlot;
            m_grants[uid] = g;
            t.ue->m_arrival += t.ue->m_prsvpSlots;
""")
open(f"{H}/ntn-nr-sidelink.cc", "w").write(s)

# ---------- example
s = open(f"{H}/orig/ntn-v2x-canyon-relay.cc").read()
s = rep(s, ' *   blindbest among neighbors with crossings apart from the source\'s, the one whose lane gives the lowest\n *             probability of being blocked toward every satellite under the street\'s height law\n',
 ' *   blindbest among neighbors with crossings apart from the source\'s, the one whose lane gives the lowest\n *             probability of being blocked toward every satellite under the street\'s height law\n'
 ' *   rand2, rand3  k distinct neighbors drawn afresh for every message\n'
 ' *   roundrobin the next neighbor of the table in cyclic order\n'
 ' *   hybrid    fresh draw among apart neighbors in the best-ranked lane (random when there is none)\n'
 ' *   anyfwd    no relay is named; every neighbor within --maxRelayM that decodes the blocked message and has an uplink forwards it\n'
 ' *   anyp      as anyfwd, but each such neighbor forwards with probability min(1, kappa / (1 + uplinks it knows within range))\n'
 ' *\n * Paper 8 revision (Nov 2026): shared per-satellite uplink (--ulKbps, --grantRtt, --ulMaxQueueS), relay link setup\n'
 ' * latency (--relaySetupMs), mobility from a trace (--mobTrace), facades from a file or with cross-street gaps\n'
 ' * (--facadeFile, --gapEveryM, --gapLenM), urban sidelink channel (--slUrban), aligned reselection (--alignResel).\n')
s = rep(s, "#include <cmath>\n#include <fstream>\n", "#include <algorithm>\n#include <cmath>\n#include <fstream>\n")
s = rep(s, "constexpr uint32_t kNoRelay = 0xFFFFFFFF;\n", "constexpr uint32_t kNoRelay = 0xFFFFFFFF;\nconstexpr uint32_t kAny = 0xFFFFFFFE;\n")
s = rep(s, "    uint32_t persist{0xFFFFFFFF}; // randpersist: the neighbor drawn once and kept while it stays usable\n",
"""    uint32_t persist{0xFFFFFFFF}; // randpersist: the neighbor drawn once and kept while it stays usable
    uint32_t rrLast{0};           // roundrobin: the neighbor named last
    uint32_t linked{0xFFFFFFFF};  // relay setup: the neighbor the source has a link with
    uint64_t blockedGen{0};       // messages generated while the uplink was down
    int dir{1};                   // +1 or -1 along the street
    std::vector<float> traceS;    // cumulative distance driven per trace sample (--mobTrace)
""")
s = rep(s, "double g_hoS = 0.0; // uplink interruption at each satellite change (s)\n",
"""double g_hoS = 0.0; // uplink interruption at each satellite change (s)
double g_ulKbps = 0;       // capacity of one satellite's uplink left for this street (kb/s); 0 = no limit
double g_grantRtt = 0;     // round trips from a scheduling request to the grant
double g_ulMaxQ = 1.0;     // largest wait for an uplink opportunity (s); later messages are dropped
double g_setupS = 0;       // latency of a link to a relay the source did not name last (s)
double g_nbAge = 1.0;      // neighbor table entries age out after this time (s)
double g_kappa = 1.0;      // anyp: expected number of forwards
double g_traceDt = 1.0;
bool g_useTrace = false;
double g_gapEvery = 0, g_gapLen = 0;
std::vector<float> g_rasterL, g_rasterR; // facades from a file, per meter along the street
std::vector<double> g_satBusy;           // per satellite: the time its uplink is next free
uint64_t g_ulDropQ = 0, g_ulDropLos = 0, g_ulSent = 0, g_ulBytes = 0;
double g_ulWaitSum = 0;
""")
s = rep(s, """    double x = std::fmod(v.x0 + v.v * Simulator::Now().GetSeconds(), g_len);
""", """    double x0 = v.x0 + v.v * Simulator::Now().GetSeconds();
    if (g_useTrace)
    {
        const double u = Simulator::Now().GetSeconds() / g_traceDt;
        const size_t k = std::min(static_cast<size_t>(u), v.traceS.size() - 2);
        const double f = std::min(1.0, u - k);
        x0 = v.x0 + v.dir * (v.traceS[k] + f * (v.traceS[k + 1] - v.traceS[k]));
    }
    double x = std::fmod(x0, g_len);
""")
s = rep(s, """    const auto& h = right ? g_hRight : g_hLeft;
    double m = std::fmod(xc, g_len);
    if (m < 0)
    {
        m += g_len;
    }
""", """    const auto& h = right ? g_hRight : g_hLeft;
    double m = std::fmod(xc, g_len);
    if (m < 0)
    {
        m += g_len;
    }
    if (!g_rasterL.empty())
    {
        const auto& ra = right ? g_rasterR : g_rasterL;
        return ra[std::min(static_cast<size_t>(m), ra.size() - 1)];
    }
    if (g_gapEvery > 0 && std::fmod(m, g_gapEvery) < g_gapLen)
    {
        return 0.0; // a cross street opens the sky
    }
""")
# uplink send path
s = rep(s, """Ptr<Packet>
MakeTx(uint32_t i)
{""", """/// One uplink message of vehicle i: a scheduling request and grant take grantRtt round trips to the satellite's
/// gateway, then the message waits its turn on the satellite's shared uplink. A vehicle that has lost its uplink by
/// the time of its turn loses the message.
void
UlFinish(uint32_t i, Ptr<Packet> p, uint16_t port, double t0)
{
    Veh& v = g_v[i];
    if (!v.up)
    {
        ++g_ulDropLos;
        return;
    }
    g_ulWaitSum += Simulator::Now().GetSeconds() - t0;
    ++g_ulSent;
    g_ulBytes += p->GetSize();
    v.sock->SendTo(p, 0, InetSocketAddress(g_srvAddr, port));
}

void
UlSend(uint32_t i, Ptr<Packet> p, uint16_t port)
{
    Veh& v = g_v[i];
    if (g_ulKbps <= 0 && g_grantRtt <= 0)
    {
        ++g_ulSent;
        g_ulBytes += p->GetSize();
        v.sock->SendTo(p, 0, InetSocketAddress(g_srvAddr, port));
        return;
    }
    const double now = Simulator::Now().GetSeconds();
    const double rtt = 2.0 * 2.0 * v.slant / kC; // service and feeder link, both ways
    double t = now + g_grantRtt * rtt;
    if (g_ulKbps > 0)
    {
        double& busy = g_satBusy[v.sat];
        t = std::max(t, busy);
        if (t - now > g_ulMaxQ)
        {
            ++g_ulDropQ;
            return;
        }
        t += (p->GetSize() + 28) * 8.0 / (g_ulKbps * 1e3);
        busy = t;
    }
    Simulator::Schedule(Seconds(t - now), &UlFinish, i, p, port, now);
}

Ptr<Packet>
MakeTx(uint32_t i)
{""")
s = rep(s, """    if (v.up)
    {
        v.sock->SendTo(body->Copy(), 0, InetSocketAddress(g_srvAddr, kPort));
        ++v.sentDirect;
    }
""", """    if (v.up)
    {
        UlSend(i, body->Copy(), kPort);
        ++v.sentDirect;
    }
    else
    {
        ++v.blockedGen;
    }
""")
s = rep(s, """                          g_policy == "randpersist";
    if (!v.up && canRelay)
    {
        relay = ChooseRelay(i);
        v.relayAsked += relay != kNoRelay;
    }
""", """                          g_policy == "randpersist" || g_policy == "rand2" || g_policy == "rand3" ||
                          g_policy == "roundrobin" || g_policy == "hybrid";
    if (!v.up && canRelay)
    {
        relay = ChooseRelay(i);
        v.relayAsked += relay != kNoRelay;
    }
    else if (!v.up && (g_policy == "anyfwd" || g_policy == "anyp"))
    {
        relay = kAny;
        ++v.relayAsked;
    }
""")
s = rep(s, """    tag.SetSeq(relay);
""", """    tag.SetSeq(relay); // kNoRelay, kAny, one id, or up to three ids: id1 + 1 + ((id2 + 1) << 10) + ((id3 + 1) << 20)
""")
# relay side
s = rep(s, """    SeqTsHeader tag;
    c->RemoveHeader(tag);
    if (tag.GetSeq() != self || !g_v[self].up)
    {
        return;
    }
    // forward the source's message unchanged; its header keeps the source's generation time
    g_v[self].sock->SendTo(c, 0, InetSocketAddress(g_srvAddr, kPort + 1));
    ++g_v[self].relayed;
}
""", """    const double fromX = ((info.GetSeq() >> 3) & 0x3FFF) / 10.0;
    const double fromY = lanesOf(info.GetSeq());
    SeqTsHeader tag;
    c->RemoveHeader(tag);
    const uint32_t tg = tag.GetSeq();
    bool named = false;
    if (tg == kAny)
    {
        const double dme = std::hypot(fromX - PosX(g_v[self]), fromY - g_v[self].y);
        if (dme > g_maxRelay)
        {
            return;
        }
        named = true;
        if (g_policy == "anyp")
        {
            const double now = Simulator::Now().GetSeconds();
            uint32_t ups = 0;
            for (const auto& [j, n] : g_v[self].nb)
            {
                ups += n.up && now - n.t <= g_nbAge && std::abs(n.x - PosX(g_v[self])) <= g_maxRelay;
            }
            static Ptr<UniformRandomVariable> up = CreateObject<UniformRandomVariable>();
            named = up->GetValue() < std::min(1.0, g_kappa / (1.0 + ups));
        }
    }
    else if (tg != kNoRelay)
    {
        for (int k = 0; k < 3; ++k)
        {
            const uint32_t id = (tg >> (10 * k)) & 0x3FF;
            named |= id != 0 && id - 1 == self;
        }
    }
    if (!named || !g_v[self].up)
    {
        return;
    }
    // forward the source's message unchanged; its header keeps the source's generation time
    if (g_setupS > 0 && tg != kAny && g_v[from].linked != self)
    {
        g_v[from].linked = self;
        Simulator::Schedule(Seconds(g_setupS), [self, c]() {
            if (g_v[self].up)
            {
                UlSend(self, c, kPort + 1);
                ++g_v[self].relayed;
            }
        });
        return;
    }
    UlSend(self, c, kPort + 1);
    ++g_v[self].relayed;
}
""")
s = rep(s, """void
SlRx(uint32_t self, uint32_t from, Ptr<Packet> p)
{""", """double
lanesOf(uint32_t s)
{
    static const double lanes[] = {4.0, 7.5, 12.5, 16.0};
    return lanes[(s >> 1) & 3u] * g_width / 20.0;
}

void
SlRx(uint32_t self, uint32_t from, Ptr<Packet> p)
{""")
# relay choice
s = rep(s, """    std::vector<uint32_t> cand;
    for (const auto& [j, n] : g_v[src].nb)
    {
        if (now - n.t > 1.0)
""", """    std::vector<uint32_t> cand;
    std::vector<std::pair<uint32_t, double>> apartQ; // hybrid: apart neighbors and the score of their lane
    const bool needApart = g_policy == "blindfar" || g_policy == "blindbest" || g_policy == "hybrid";
    for (const auto& [j, n] : g_v[src].nb)
    {
        if (now - n.t > g_nbAge)
""")
s = rep(s, """        bool apart = true;
        for (const auto& [e1, a1] : g_visEA)
        {
            for (const auto& [e2, a2] : g_visEA)
            {
""", """        bool apart = needApart;
        for (const auto& [e1, a1] : g_visEA)
        {
            if (!needApart)
            {
                break;
            }
            for (const auto& [e2, a2] : g_visEA)
            {
""")
s = rep(s, """            if (q < bestQ - 1e-9 || (std::abs(q - bestQ) <= 1e-9 && d < bestQD))
""", """            apartQ.emplace_back(j, q);
            if (q < bestQ - 1e-9 || (std::abs(q - bestQ) <= 1e-9 && d < bestQD))
""")
s = rep(s, """    if (g_policy == "blindfar")
    {
        if (pickFar != kNoRelay)
""", """    if ((g_policy == "rand2" || g_policy == "rand3") && !cand.empty())
    {
        static Ptr<UniformRandomVariable> uk = CreateObject<UniformRandomVariable>();
        const uint32_t k = g_policy == "rand2" ? 2 : 3;
        uint32_t tagv = 0;
        for (uint32_t m = 0; m < k && !cand.empty(); ++m)
        {
            const uint32_t idx = uk->GetInteger(0, cand.size() - 1);
            tagv |= (cand[idx] + 1) << (10 * m);
            cand.erase(cand.begin() + idx);
        }
        return tagv;
    }
    if (g_policy == "roundrobin" && !cand.empty())
    {
        auto& v = g_v[src];
        auto it = std::upper_bound(cand.begin(), cand.end(), v.rrLast); // the table is ordered by id
        v.rrLast = it == cand.end() ? cand.front() : *it;
        return v.rrLast + 1;
    }
    if (g_policy == "hybrid" && !cand.empty())
    {
        static Ptr<UniformRandomVariable> uh = CreateObject<UniformRandomVariable>();
        std::vector<uint32_t> best;
        for (const auto& [j, q] : apartQ)
        {
            if (q <= bestQ + 1e-9)
            {
                best.push_back(j);
            }
        }
        if (best.empty())
        {
            ++g_farFallback;
            return cand[uh->GetInteger(0, cand.size() - 1)] + 1;
        }
        return best[uh->GetInteger(0, best.size() - 1)] + 1;
    }
    if (g_policy == "blindfar")
    {
        if (pickFar != kNoRelay)
""")
# single-id returns become id + 1, keeping kNoRelay as is
s = rep(s, "            return pickFar;\n", "            return pickFar + 1;\n")
s = rep(s, "        return v.persist;\n", "        return v.persist + 1;\n")
s = rep(s, "            return pickBest;\n", "            return pickBest + 1;\n")
s = rep(s, "        return pickLb; // ablation", "        return pickLb + 1; // ablation")
s = rep(s, "        return cand[u->GetInteger(0, cand.size() - 1)]; // ablation", "        return cand[u->GetInteger(0, cand.size() - 1)] + 1; // ablation")
s = rep(s, "    return pick;\n}\n\nPtr<Packet>" if False else "    return pick;\n}\n", "    return pick == kNoRelay ? kNoRelay : pick + 1;\n}\n")
# main
s = rep(s, """    cmd.AddValue("lat", "street latitude", g_lat);
""", """    cmd.AddValue("ulKbps", "capacity of a satellite's uplink left for this street (kb/s); 0 = no limit", g_ulKbps);
    cmd.AddValue("grantRtt", "round trips from scheduling request to grant", g_grantRtt);
    cmd.AddValue("ulMaxQueueS", "largest wait for an uplink opportunity (s)", g_ulMaxQ);
    double setupMs = 0;
    cmd.AddValue("relaySetupMs", "latency of a link to a relay the source did not name last (ms)", setupMs);
    cmd.AddValue("nbAgeS", "neighbor entries age out after this time (s)", g_nbAge);
    cmd.AddValue("kappa", "anyp: expected number of forwards", g_kappa);
    std::string mobTrace, facadeFile;
    cmd.AddValue("mobTrace", "CSV: first line dt, then one line per vehicle: x0,dir,s(0),s(dt),...", mobTrace);
    cmd.AddValue("facadeFile", "CSV: side (0 left, 1 right),x0,x1,height", facadeFile);
    cmd.AddValue("gapEveryM", "a cross street every this many meters (0 = none)", g_gapEvery);
    cmd.AddValue("gapLenM", "width of a cross street (m)", g_gapLen);
    bool slUrban = false, alignResel = false;
    double slNlosP = 0.0;
    cmd.AddValue("slNlosP", "stress: share of vehicle pairs given the NLOS law", slNlosP);
    cmd.AddValue("slUrban", "TR 37.885 urban sidelink channel", slUrban);
    cmd.AddValue("alignResel", "message arrivals on a fixed period with a random phase", alignResel);
    cmd.AddValue("lat", "street latitude", g_lat);
""")
s = rep(s, """    g_hSig = hSigma;
""", """    g_hSig = hSigma;
    g_setupS = setupMs / 1000.0;
    if (!facadeFile.empty())
    {
        std::ifstream f(facadeFile);
        NS_ABORT_MSG_UNLESS(f.good(), "cannot open " << facadeFile);
        const size_t nm = static_cast<size_t>(std::ceil(g_len));
        g_rasterL.assign(nm, 0.0f);
        g_rasterR.assign(nm, 0.0f);
        std::string line;
        while (std::getline(f, line))
        {
            int side;
            double a, b, hh;
            if (std::sscanf(line.c_str(), "%d,%lf,%lf,%lf", &side, &a, &b, &hh) != 4)
            {
                continue;
            }
            auto& ra = side ? g_rasterR : g_rasterL;
            for (size_t m = static_cast<size_t>(std::max(0.0, a)); m < nm && m < b; ++m)
            {
                ra[m] = static_cast<float>(hh);
            }
        }
    }
""")
s = rep(s, """            g_sats.push_back(m);
        }
    }
""", """            g_sats.push_back(m);
        }
    }
    g_satBusy.assign(g_sats.size(), 0.0);
""")
s = rep(s, """    Ptr<NtnSlChannel> ch = CreateObject<NtnSlChannel>();
""", """    Ptr<NtnSlChannel> ch = CreateObject<NtnSlChannel>();
    ch->SetUrbanChannel(slUrban, slNlosP);
""")
s = rep(s, """        v.mac->SetReservationPeriod(periodMs);
""", """        v.mac->SetReservationPeriod(periodMs);
        v.mac->SetAlignReselection(alignResel);
""")
s = rep(s, """    for (uint16_t via = 0; via < 2; ++via)
    {
""", """    if (!mobTrace.empty())
    {
        std::ifstream f(mobTrace);
        NS_ABORT_MSG_UNLESS(f.good(), "cannot open " << mobTrace);
        std::string line;
        std::getline(f, line);
        g_traceDt = std::stod(line);
        for (uint32_t i = 0; i < nVeh; ++i)
        {
            NS_ABORT_MSG_UNLESS(std::getline(f, line), "trace has too few vehicles");
            std::stringstream ls(line);
            std::string tok;
            std::getline(ls, tok, ',');
            g_v[i].x0 = std::stod(tok);
            std::getline(ls, tok, ',');
            g_v[i].dir = std::stoi(tok);
            while (std::getline(ls, tok, ','))
            {
                g_v[i].traceS.push_back(std::stof(tok));
            }
        }
        g_useTrace = true;
    }
    for (uint16_t via = 0; via < 2; ++via)
    {
""")
s = rep(s, """    uint64_t gen = 0, direct = 0, relayed = 0, asked = 0;
""", """    uint64_t gen = 0, direct = 0, relayed = 0, asked = 0, blockedGen = 0;
    for (const auto& v : g_v)
    {
        blockedGen += v.blockedGen;
    }
""")
s = rep(s, """              << " relayed=" << relayed;
""", """              << " relayed=" << relayed << " blocked_gen=" << blockedGen << " ul_sent=" << g_ulSent
              << " ul_drop_q=" << g_ulDropQ << " ul_drop_los=" << g_ulDropLos << " ul_bytes=" << g_ulBytes
              << " ul_wait_ms=" << (g_ulSent ? 1e3 * g_ulWaitSum / g_ulSent : 0.0);
""")
s = rep(s, '#include "ns3/ntn-nr-sidelink.h"', '#include "ntn-nr-sidelink.h"')
open(f"{H}/ntn-v2x-canyon-relay.cc", "w").write(s)
