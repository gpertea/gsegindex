#!/bin/env bash
tparms=' <<<<\nWallTime: %E\nUserTime: %U\n Max Mem: %M\n%%CPU : %P'
/usr/bin/time "-f>>>> ailist$tparms" ./ailist d4.bed d0.bed > ailist.out
/usr/bin/time "-f>>>> gailist$tparms" ./gailist d4.bed d0.bed > gailist.out
diff gailist.out ailist.out | head -20
