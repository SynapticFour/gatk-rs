import com.intel.gkl.pairhmm.IntelPairHmm;
import java.io.BufferedReader;
import java.io.FileReader;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.List;
import org.broadinstitute.gatk.nativebindings.pairhmm.HaplotypeDataHolder;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.gatk.nativebindings.pairhmm.ReadDataHolder;

/** Stock GATK 4.4.0.0 libgkl_pairhmm.so likelihoods for the frozen 6R.284 inputs. */
public final class LiveGkl284 {
    private static byte[] parseBytes(String csv) {
        if (csv.isEmpty()) {
            return new byte[0];
        }
        String[] parts = csv.split(",");
        byte[] out = new byte[parts.length];
        for (int i = 0; i < parts.length; i++) {
            out[i] = (byte) Integer.parseInt(parts[i]);
        }
        return out;
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 1) {
            throw new IllegalArgumentException("usage: LiveGkl284 inputs.tsv");
        }
        List<String> hapBases = new ArrayList<String>();
        List<Integer> readIndex = new ArrayList<Integer>();
        List<Integer> skipped = new ArrayList<Integer>();
        List<ReadDataHolder> reads = new ArrayList<ReadDataHolder>();
        try (BufferedReader br = new BufferedReader(new FileReader(args[0]))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '#') {
                    continue;
                }
                String[] f = line.split("\t", -1);
                if (f[0].equals("HAP")) {
                    hapBases.add(f[3]);
                } else if (f[0].equals("READ")) {
                    ReadDataHolder r = new ReadDataHolder();
                    r.readBases = f[6].getBytes("US-ASCII");
                    r.readQuals = parseBytes(f[7]);
                    r.insertionGOP = parseBytes(f[8]);
                    r.deletionGOP = parseBytes(f[9]);
                    r.overallGCP = parseBytes(f[10]);
                    reads.add(r);
                    readIndex.add(Integer.valueOf(f[1]));
                    skipped.add(Integer.valueOf(f[4]));
                }
            }
        }
        HaplotypeDataHolder[] haps = new HaplotypeDataHolder[hapBases.size()];
        for (int i = 0; i < haps.length; i++) {
            haps[i] = new HaplotypeDataHolder();
            haps[i].haplotypeBases = hapBases.get(i).getBytes("US-ASCII");
        }
        ReadDataHolder[] readArr = reads.toArray(new ReadDataHolder[0]);
        double[] like = new double[readArr.length * haps.length];
        IntelPairHmm hmm = new IntelPairHmm();
        if (!hmm.load(null)) {
            throw new IllegalStateException("stock libgkl_pairhmm.so did not load");
        }
        PairHMMNativeArguments nativeArgs = new PairHMMNativeArguments();
        nativeArgs.useDoublePrecision = false;
        nativeArgs.maxNumberOfThreads = 1;
        hmm.initialize(nativeArgs);
        hmm.computeLikelihoods(readArr, haps, like);
        hmm.done();
        int n = 0;
        for (int r = 0; r < readArr.length; r++) {
            for (int h = 0; h < haps.length; h++) {
                double v = like[n++];
                long bits = Double.doubleToRawLongBits(v);
                System.out.printf(
                        "STOCK\t%d\t%d\t%d\t0x%016x\t%.17g%n",
                        readIndex.get(r).intValue(),
                        h,
                        skipped.get(r).intValue(),
                        bits,
                        v);
            }
        }
    }
}
