#!/bin/bash
set -euo pipefail

# ---------- config ----------
R_SCRIPT="code_1_distance.R"
N_TASKS=10
TASK_NAME="Heterogeneity-Digits-PairwiseDistances"
# --------------------------

# Avoid oversubscription (1 thread per R process)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1

# Email config
RECIPIENT="kisung.you@outlook.com"
FROM_NAME="Kisung-VC"
FROM_ADDR="ksyou496@gmail.com"

start_ts=$(date +%s)
host=$(hostname -f 2>/dev/null || hostname)

# cores/jobs (manual)
CORES=12
JOBS=$((CORES - 2))
JOBS=$((JOBS > 0 ? JOBS : 1))
JOBS=$((JOBS > N_TASKS ? N_TASKS : JOBS))

echo "Running with $JOBS parallel jobs (cores=$CORES, tasks=$N_TASKS)"

# ---- PARALLEL EXECUTION (LIVE OUTPUT) ----
seq 1 "$N_TASKS" | parallel --jobs "$JOBS" --halt soon,fail=1 '
  echo "[START $(date +%H:%M:%S)] task={}"
  Rscript "'"$R_SCRIPT"'" {}
  rc=$?
  echo "[END   $(date +%H:%M:%S)] task={} rc=$rc"
  exit $rc
'
parallel_exit=$?

# -----------------------------------------

end_ts=$(date +%s)
elapsed=$((end_ts - start_ts))

h=$((elapsed/3600))
m=$(((elapsed%3600)/60))
s=$((elapsed%60))
elapsed_hms=$(printf "%02d:%02d:%02d" "$h" "$m" "$s")

status="OK"
if (( parallel_exit != 0 )); then
  status="FAILED (exit $parallel_exit)"
fi

subject="${TASK_NAME} - R batch ${status} on ${host} — ${N_TASKS} tasks, ${JOBS} jobs"
body=$(cat <<EOF
Host      : ${host}
Script    : ${R_SCRIPT}
Tasks     : ${N_TASKS}
Parallel  : ${JOBS}
Started   : $(date -d "@${start_ts}" "+%Y-%m-%d %H:%M:%S %Z")
Finished  : $(date -d "@${end_ts}"   "+%Y-%m-%d %H:%M:%S %Z")
Elapsed   : ${elapsed_hms}
Result    : ${status}

This message was sent automatically after GNU Parallel finished.
EOF
)

printf "%s\n" "$body" \
| mail -aFrom:"\"${FROM_NAME}\" <${FROM_ADDR}>" \
       -s "$subject" "$RECIPIENT"

echo
echo "$subject"
echo "$body"
exit "$parallel_exit"
