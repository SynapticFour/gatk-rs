#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
matrix="${repo_root}/parity/fixtures/p6_pairhmm_live_matrix.tsv"
repo_container="/work"
matrix_container="${repo_container}/parity/fixtures/p6_pairhmm_live_matrix.tsv"
mkdir -p "${report_dir}"

gatk_image="${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}"
gatk_platform="${GATK_DOCKER_PLATFORM:-linux/amd64}"
warn_threshold="${P6_DRIFT_WARN_THRESHOLD:-1.0}"
fail_threshold="${P6_DRIFT_FAIL_THRESHOLD:-5.0}"
rust_gap_open="${P6_GAP_OPEN_PROB:-0.005}"
rust_gap_extend="${P6_GAP_EXTEND_PROB:-0.1}"
rust_ins_emit="${P6_INSERTION_EMISSION_PROB:-0.25}"

java_out="${report_dir}/p6_pairhmm_live.java.tsv"
rust_out="${report_dir}/p6_pairhmm_live.rust.tsv"
java_out_container="${repo_container}/parity/reports/p6_pairhmm_live.java.tsv"
summary_json="${report_dir}/p6_pairhmm_live_drift_summary.json"
summary_md="${report_dir}/p6_pairhmm_live_drift_summary.md"

echo "[p6-live-drift] matrix=${matrix}"
echo "[p6-live-drift] java=${gatk_image} (${gatk_platform})"
echo "[p6-live-drift] rust params: gap_open=${rust_gap_open} gap_extend=${rust_gap_extend} ins_emit=${rust_ins_emit}"

docker run --rm --platform "${gatk_platform}" "${gatk_image}" gatk --version >/dev/null

docker run --rm --platform "${gatk_platform}" \
  -v "${repo_root}:/work" -w /work \
  "${gatk_image}" \
  bash -lc "cat > /tmp/P6PairHmmDump.java <<'EOF'
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import org.broadinstitute.hellbender.utils.pairhmm.Log10PairHMM;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;

public class P6PairHmmDump {
  static byte[] parseQuals(String raw) {
    String[] parts = raw.split(\",\");
    byte[] out = new byte[parts.length];
    for (int i = 0; i < parts.length; i++) out[i] = (byte) Integer.parseInt(parts[i]);
    return out;
  }

  static byte[] fill(int n, int value) {
    byte[] out = new byte[n];
    Arrays.fill(out, (byte) value);
    return out;
  }

  public static void main(String[] args) throws Exception {
    if (args.length != 2) throw new IllegalArgumentException(\"usage: <input.tsv> <output.tsv>\");
    List<String> lines = Files.readAllLines(Paths.get(args[0]));
    Method method = PairHMM.class.getDeclaredMethod(
        \"computeReadLikelihoodGivenHaplotypeLog10\",
        byte[].class, byte[].class, byte[].class, byte[].class, byte[].class, byte[].class,
        boolean.class, byte[].class);
    method.setAccessible(true);
    Log10PairHMM hmm = new Log10PairHMM(true);
    StringBuilder out = new StringBuilder();
    out.append(\"# case_id\\tlog10_likelihood\\n\");
    for (String line : lines) {
      if (line.isBlank() || line.startsWith(\"#\")) continue;
      String[] c = line.split(\"\\t\");
      if (c.length < 5) continue;
      String caseId = c[0];
      byte[] read = c[1].getBytes();
      byte[] baseQuals = parseQuals(c[2]);
      byte[] hap = c[4].getBytes();
      byte[] insQuals = fill(read.length, 45);
      byte[] delQuals = fill(read.length, 45);
      byte[] gcps = fill(read.length, 10);
      hmm.initialize(read.length, hap.length);
      double ll = (double) method.invoke(hmm, hap, read, baseQuals, insQuals, delQuals, gcps, true, null);
      out.append(caseId).append('\\t').append(ll).append('\\n');
    }
    hmm.close();
    Files.writeString(Paths.get(args[1]), out.toString());
  }
}
EOF
javac -cp /gatk/gatk-package-4.4.0.0-local.jar /tmp/P6PairHmmDump.java
java -cp /gatk/gatk-package-4.4.0.0-local.jar:/tmp P6PairHmmDump '${matrix_container}' '${java_out_container}'"

echo "[p6-live-drift] generating Rust likelihood vector"
CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}" RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}" \
  P6_GAP_OPEN_PROB="${rust_gap_open}" \
  P6_GAP_EXTEND_PROB="${rust_gap_extend}" \
  P6_INSERTION_EMISSION_PROB="${rust_ins_emit}" \
  cargo run -p gatk-haplotypecaller --example p6_pairhmm_dump --locked -- \
  "${matrix}" "${rust_out}" >/dev/null

echo "[p6-live-drift] comparing drift"
python3 "${repo_root}/scripts/parity/compare_p6_pairhmm_drift.py" \
  --matrix "${matrix}" \
  --java "${java_out}" \
  --rust "${rust_out}" \
  --json-out "${summary_json}" \
  --md-out "${summary_md}" \
  --warn-threshold "${warn_threshold}" \
  --fail-threshold "${fail_threshold}" \
  --gap-open "${rust_gap_open}" \
  --gap-extend "${rust_gap_extend}" \
  --ins-emission "${rust_ins_emit}"

echo "[p6-live-drift] wrote ${summary_json}"
echo "[p6-live-drift] wrote ${summary_md}"
