package FunctionV;

import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.List;

public class CountSite {

    /**
     * 对VCF文件进行变异统计.
     *
     * @param infileDirS
     */
    public static void countSitesinFastCallformat_fromTxt(String infileDirS) {

        List<File> fsList = AoFile.getFileListInDir(infileDirS);

        System.out.println("Chr\tN_RawSNPs\tN_BiallelicSNPs\tN_TriallelicSNPs\tIndels\tInsertions\tDeletions");
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
            System.out.println(String.valueOf(chrint) + "\t" + String.valueOf(cntSNP) + "\t" + String.valueOf(cntBi) + "\t" + String.valueOf(cntTri) + "\t" + String.valueOf(cntIndel) + "\t" + String.valueOf(cntI) + "\t" + String.valueOf(cntD));

        }
//        fsList.parallelStream().forEach(f -> {
//            try {
//                String infileS = f.getAbsolutePath();
//                BufferedReader br = AoFile.readFile(infileS);
//                String chr = f.getName().substring(3, 6); //提取染色体号 001
//                int chrint = Integer.parseInt(chr); //将染色体号转化为数字
//                int cntSNP = 0;
//                int cntBi = 0;
//                int cntTri = 0;
//                int cntIndel = 0;
//                int cntI = 0;
//                int cntD = 0;
//                String temp = null;
//                String header = br.readLine();
//                while ((temp = br.readLine()) != null) { //是否含有D I
//                    String alt = PStringUtils.fastSplit(temp).get(3);
//                    /**
//                     * 含有逗号，即有2个alt.
//                     * case 1: A,T
//                     * case 2: D,I
//                     * case 3: D,G
//                     * case 4: I,T
//                     */
//                    if (!(alt.length() == 1)) { //2个alt的情况;若该位点含有D或I ，那么就属于Indel，如果没有D 或者I，那么就属于SNP
//                        boolean ifD = false;
//                        if (!alt.contains("D") && (!alt.contains("I"))) {
//                            cntTri++;
//                            cntSNP++;
//                        }
//                        if (alt.contains("D")) {
//                            cntD++;
//                            cntIndel++;
//                            ifD = true;
//                        }
//                        if (alt.contains("I")) {
//                            cntI++;
//                            if (ifD == false) { //针对 case 2的情况，即含有D又含有I， 这时在上面的D中已经加过 cntIndel了，所以不用再加了
//                                cntIndel++;
//                            }
//                        }
//
//                    } else if (alt.length() == 1) { //1个alt的情况;
//                        if (!alt.equals("D") && (!alt.equals("I"))) {
//                            cntBi++;
//                            cntSNP++;
//                        }
//                        if (alt.equals("D")) {
//                            cntD++;
//                            cntIndel++;
//                        }
//                        if (alt.equals("I")) {
//                            cntI++;
//                            cntIndel++;
//                        }
//                    }
//                }
//                br.close();
//                System.out.println(String.valueOf(chrint) + "\t" + String.valueOf(cntSNP) + "\t" + String.valueOf(cntBi) + "\t" + String.valueOf(cntTri) + "\t" + String.valueOf(cntIndel) + "\t" + String.valueOf(cntI) + "\t" + String.valueOf(cntD) + "\t" + repeatDI);
//            } catch (Exception e) {
//                e.printStackTrace();
//                System.exit(1);
//            }
//        });
    }
}
