#!/bin/bash
# after (10, 13): (10, 15) in the Section 5.1 setting, then the zymo-preset runs
until [ -f ~/zymo/rb/run3b.log ] && grep -q '^== done' ~/zymo/rb/run3b.log && ! systemctl --user is-active -q zy-rb-run5 && ! systemctl --user is-active -q zy-rb-run4; do sleep 30; done
echo "== (10, 15) $(date +%T)"; bash ~/zymo/rb/run3c.sh > ~/zymo/rb/run3c.log 2>&1; echo "exit $?"
echo "== zymo preset runs $(date +%T)"; bash ~/zymo/kebab2/fixedL.sh > ~/zymo/kebab2/fixedL.log 2>&1; echo "exit $?"
echo "== chain done $(date +%T)"
