package FunctionV;

import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CountSite {

    /**
     * 对VCF文件进行变异统计.
     *
     * @param infileDirS
     */
    public static void countSitesinFastCallformat_fromTxt(String infileDirS, String outfile) {

        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        BufferedWriter bw = null;
        bw = IOUtils.getTextWriter(outfile);
        try {
            bw.write("Chr\tN_RawSNPs\tN_BiallelicSNPs\tN_TriallelicSNPs\tIndels\tInsertions\tDeletions\n");
        } catch (IOException e) {
            e.printStackTrace();
        }
//        System.out.println("Chr\tN_RawSNPs\tN_BiallelicSNPs\tN_TriallelicSNPs\tIndels\tInsertions\tDeletions");
        for (File f : fsList){
            String infileS = f.getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);
            String chr = f.getName().substring(3, 6); //提取染色体号 001
            int chrint = Integer.parseInt(chr); //将染色体号转化为数字
            int cntSNP = 0;
            int cntBi = 0;
            int cntTri = 0;
            int cntIndel = 0;
            int cntI = 0;
            int cntD = 0;
            String temp = null;
            try {
                String header = br.readLine();
            } catch (IOException e) {
                e.printStackTrace();
            }
            while (true) {
                try {
                    if (!((temp = br.readLine()) != null)) break;
                    if (temp.startsWith("#")) {
                        continue;
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                } //是否含有D I
                String alt = PStringUtils.fastSplit(temp).get(4);
                /**
                 * 含有逗号，即有2个alt.
                 * case 1: A,T
                 * case 2: D,I
                 * case 3: D,G
                 * case 4: I,T
                 */
                if (!(alt.length() == 1)) { //2个alt的情况;若该位点含有D或I ，那么就属于Indel，如果没有D 或者I，那么就属于SNP
                    boolean ifD = false;
                    if (!alt.contains("D") && (!alt.contains("I"))) {
                        cntTri++;
                        cntSNP++;
                    }
                    if (alt.contains("D")) {
                        cntD++;
                        cntIndel++;
                        ifD = true;
                    }
                    if (alt.contains("I")) {
                        cntI++;
                        if (ifD == false) { //针对 case 2的情况，即含有D又含有I， 这时在上面的D中已经加过 cntIndel了，所以不用再加了
                            cntIndel++;
                        }
                    }

                } else if (alt.length() == 1) { //1个alt的情况;
                    if (!alt.equals("D") && (!alt.equals("I"))) {
                        cntBi++;
                        cntSNP++;
                    }
                    if (alt.equals("D")) {
                        cntD++;
                        cntIndel++;
                    }
                    if (alt.equals("I")) {
                        cntI++;
                        cntIndel++;
                    }
                }
            }
            try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
            try {
                bw.write(String.valueOf(chrint) + "\t" + String.valueOf(cntSNP) + "\t" + String.valueOf(cntBi) + "\t" + String.valueOf(cntTri) + "\t" + String.valueOf(cntIndel) + "\t" + String.valueOf(cntI) + "\t" + String.valueOf(cntD)+"\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
//            System.out.println(String.valueOf(chrint) + "\t" + String.valueOf(cntSNP) + "\t" + String.valueOf(cntBi) + "\t" + String.valueOf(cntTri) + "\t" + String.valueOf(cntIndel) + "\t" + String.valueOf(cntI) + "\t" + String.valueOf(cntD));

        }
        try {
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    /**
     * 将raw vcf文件中的snp和indels分开输出到两个文件中。
     */
    public static void splitVcf(String infileDirS) {
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        for (File f : fsList){
            String temp = null;
            BufferedWriter bwID = null;
            BufferedWriter bwSNP = null;
            String infileS = f.getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);

            String outfileID = infileS+".indel.gz";
            String outfileSNP = infileS+".snp.gz";
            bwID = IOUtils.getTextGzipWriter(outfileID);
            bwSNP = IOUtils.getTextGzipWriter(outfileSNP);
            while (true) {
                try {
                    if (!((temp = br.readLine()) != null)) break;
                    if (temp.startsWith("#")) {
                        bwID.write(temp+"\n");
                        bwSNP.write(temp+"\n");
                        continue;
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
                String alt = PStringUtils.fastSplit(temp).get(4);
                if (alt.contains("D") || (alt.contains("I"))) {
                    try {
                        bwID.write(temp+"\n");
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                } else {
                    try {
                        bwSNP.write(temp+"\n");
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
            try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
            try {
                bwID.close();
                bwSNP.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
    /**
     * 结果只有一列 depth平均深度
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public static void calVcfAverageDepth(String infileDirS, String outfileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = null;
            BufferedReader br = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".depth.txt.gz")).getAbsolutePath();
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".depth.txt.gz")).getAbsolutePath();
            }
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                while ((temp = br.readLine()).startsWith("##")) {
                }
                List<String> linetaxa = PStringUtils.fastSplit(temp, "\t");
                bw.write("AverageDepth");
                bw.newLine();
                String[] taxa = new String[linetaxa.size() - 9];
                cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++; // 对snp开始计数
                    if (cnt % 1000000 == 0) {
                        System.out.println(String.valueOf(cnt) + " lines");
                    }

                    TDoubleArrayList depthList = new TDoubleArrayList();
                    TDoubleArrayList PValueList = new TDoubleArrayList();
                    List<String> l = new ArrayList<>();
                    l = PStringUtils.fastSplit(temp, "\t");
                    String chr = l.get(0);
                    String pos = l.get(1);
                    for (int i = 0; i < taxa.length; i++) {
                        String genoS = l.get(i + 9);
                        if (genoS.startsWith(".")) {
                            depthList.add(0);
                            continue;
                        }
                        List<String> ll = PStringUtils.fastSplit(genoS, ":");
                        List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                        int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1));
                        depthList.add(depth);
                    }
                    double[] dep = depthList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(dep);
                    double relativeMean = d.getMean();
                    double sd = d.getStandardDeviation();
                    StringBuilder sb = new StringBuilder();
                    sb.append(String.format("%.4f", relativeMean));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
                System.out.println();
                System.out.println(f.getName() + "\t" + cnt + " sites is completed");
            } catch (Exception e) {
                System.out.println(temp);
                e.printStackTrace();
                System.exit(1);
            }

        });
    }
}
