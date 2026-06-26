# GRETAAnalysis (TrackBackToMe)

Same codebase as `/Users/heather/GRETAAnalysis`, with TrackBackToMe-specific settings (`EB_DIFF_TIME=210`, crystal remap in `GRETA.cpp`).

See `/Users/heather/GRETAAnalysis/README.md` for full build and workflow docs. Quick reference:

```bash
scons readGreta trackGreta

# Pass 1
./readGreta -f <data> -rootFile out.root

# Pass 2
./trackGreta -i out.root -o tracked.root -trackingChat tracking/greta_readGreta.chat

# Parallel pass 2
./scripts/trackGretaParallel.sh -j 8 -i out.root -o tracked.root \
    -trackingChat tracking/greta_readGreta.chat
```
