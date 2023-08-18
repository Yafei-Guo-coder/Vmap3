package MainV;


import FunctionV.CountSite;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import pgl.infra.utils.IOUtils;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import static FunctionV.CountSite.splitVcf;

public class Start {
//    public static void test() throws IOException {
    public static void main(String args[]) throws IOException {

//        String a = "A";
//        System.out.println(a.length());
        splitVcf(args[0]);
        //splitVcf("/data2/yafei/004_Vmap3/VCF/Raw_VCF/AA_vcf/chr019");
//        String inputDirS = "/Users/guoyafei/Documents/02_VmapIII/03_Fastcall2/测试数据/fastCall/ing/D_0026/43_1_452529.ing.gz";
////        String inputDirS = input;
////        String outputDirS = "/Users/guoyafei/Desktop/Out/43_1_452529.ing.gz";
//        DataInputStream dis1 = IOUtils.getBinaryGzipReader(inputDirS);
//        IntArrayList codedAlleleInfo = null;
//        IntArrayList codedAlleleInfo2 = null;
//        String taxonName = dis1.readUTF();
//        short chrom = dis1.readShort();
//        int binStart = dis1.readInt();
//        int binEnd = dis1.readInt();
//
//        codedAlleleInfo = new IntArrayList();
//        int currentRecord = 0;
//        while ((currentRecord = dis1.readInt()) != Integer.MIN_VALUE) {
//            int v = currentRecord >> 9;
//            codedAlleleInfo.add(v);
//            System.out.println(v+"\n");
//            break;
//        }
//        System.out.println(codedAlleleInfo.size());
//        DataInputStream dis2 = IOUtils.getBinaryGzipReader(outputDirS);
//        String taxonName2 = dis2.readUTF();
//        short chrom2 = dis2.readShort();
//        int binStart2 = dis2.readInt();
//        int binEnd2 = dis2.readInt();

//        codedAlleleInfo2 = new IntArrayList();
//        int currentRecord2 = 0;
//        while ((currentRecord2 = dis2.readInt()) != Integer.MIN_VALUE) {
//            System.out.println(codedAlleleInfo2.size());
//            codedAlleleInfo2.add(currentRecord2);
//        }
//        System.out.println(codedAlleleInfo2.size());
//        System.out.println(codedAlleleInfo.size());

//        BufferedWriter bw = IOUtils.getTextWriter(this.vLibPosFileS);
//        for (int i = 0; i < ; i++) {
//            sb.setLength(0);
//            sb.append(this.chrom).append("\t").append(vl.getPosition(i));
//            bw.write(sb.toString());
//            bw.newLine();
//        }
    }
//    public static void main(String args[]) throws IOException {
//        List<File> Dirlist = IOUtils.getFileListInDir("/Users/guoyafei/Desktop/ing");
//
////        for (int j = 0; j < Dirlist.size(); j++) {
//            List<File> Dir2list = IOUtils.getFileListInDir(String.valueOf(Dirlist.get(25)));
//                test(Dir2list.get(0).getAbsolutePath());
//            getHaplo(Dir2list.get(1).getAbsolutePath(),Dirlist.get(j).getAbsolutePath().replace("vcf","txt"));
//        }
//        CountSite.countSitesinFastCallformat_fromTxt(args[0],args[1]);
//        String[] taxaNames = null;
//        String inputDirS = "/data2/xinyue/Yafei/vmap3_E1/AABB_middle/ing";
//        String outputDirS = "/data2/xinyue/Yafei/vmap3_E1/AABB_middle/ing2";
//        File outDir = new File (outputDirS);
////        outDir.mkdir();
//        List<File> chrDirList = IOUtils.getDirListInDir(inputDirS);
//        File[] chrOutDirs = new File[chrDirList.size()];
//        for (int i = 0; i < chrDirList.size(); i++) {
//            File f = new File (outDir, chrDirList.get(i).getName());
//            f.mkdir();
//            chrOutDirs[i] = f;
//        }
//        for (int i = 0; i < chrDirList.size(); i++) {
//            List<File> taxaDirlist = IOUtils.getDirListInDir(String.valueOf(chrDirList.get(i)));
//            taxaNames = new String[taxaDirlist.size()];
//            for (int j = 0; j < taxaDirlist.size(); j++) {
//                File f = new File(chrOutDirs[i], taxaDirlist.get(j).getName());
//                f.mkdir();
//            }
//        }
//        for (int i = 0; i < chrDirList.size(); i++) {
//            List<File> taxaDirlist = IOUtils.getDirListInDir(String.valueOf(chrDirList.get(i)));
//            taxaNames = new String[taxaDirlist.size()];
//            for (int j = 0; j < taxaDirlist.size(); j++) {
//                taxaNames[j] = taxaDirlist.get(j).getName();
//            }
//            Arrays.sort(taxaNames);
//            for (int k = 0; k < taxaDirlist.size(); k++) {
////                String taxaS = new File(String.valueOf(taxaDirlist.get(k))).getAbsolutePath();
////                String outtaxaS = taxaS.replace("ing", "ing2");
////                File f = new File (outtaxaS);
////                f.mkdir();
//                List<File> ingList = Arrays.asList(taxaDirlist.get(k).listFiles());
//                for (int j = 0; j < ingList.size(); j++) {
//                    String fileS = new File(String.valueOf(ingList.get(j))).getAbsolutePath();
////                        System.out.println(fileS);
////                    IndividualGenotype ing = new IndividualGenotype(fileS);
//                    DataInputStream dis = IOUtils.getBinaryGzipReader(fileS);
//                    String taxonName = dis.readUTF();
//                    short chrom = dis.readShort();
//                    int binStart = dis.readInt();
//                    int binEnd = dis.readInt();
//                    IntArrayList codedAlleleInfo = null;
//                    codedAlleleInfo = new IntArrayList();
//                    int currentRecord = 0;
//                    while ((currentRecord = dis.readInt()) != Integer.MIN_VALUE) {
//                        codedAlleleInfo.add(currentRecord);
//                    }
//                    for (int m = 0; m < codedAlleleInfo.size(); m++) {
//                        int v = codedAlleleInfo.getInt(m) & 192;
//                        if (v != 192) continue;
////                            System.out.println(Integer.toBinaryString(codedAlleleInfo.getInt(m)));
////                            System.out.println(ing.getAlleleBase(m));
////                            v = codedAlleleInfo.getInt(m) & 2147483392;
////                            v = v + 191;
//                        codedAlleleInfo.set(m,codedAlleleInfo.getInt(m)-1);
////                            ing.codedAlleleInfo.set(m,codedAlleleInfo.getInt(m)-1);
//
////                            System.out.println(Integer.toBinaryString(codedAlleleInfo.getInt(m)));
////                            System.out.println(ing.getAlleleBase(m));
////                            System.out.println("------------------");
//                    }
//                    String sb = new File(String.valueOf(ingList.get(j))).getAbsolutePath();
//                    String sb2 = sb.replaceFirst("ing", "ing2");
//                    DataOutputStream dos = IOUtils.getBinaryGzipWriter(sb2);
//                    try {
//                        dos.writeUTF(taxonName);
//                        dos.writeShort((short)chrom);
//                        dos.writeInt(binStart);
//                        dos.writeInt(binEnd);
//                    }
//                    catch (Exception e) {
//                        e.printStackTrace();
//                        System.exit(1);
//                    }
//                    for (int m = 0; m < codedAlleleInfo.size(); m++) {
//                        dos.writeInt(codedAlleleInfo.getInt(m));
//                    }
//                    dos.writeInt(Integer.MIN_VALUE);
//                    dos.flush();
//                    dos.close();
//                    dis.close();
//                }
////
//            }
//        splitVcf(args[0]);
//            splitVcf("/data2/yafei/004_Vmap3/VCF/Raw_VCF/AA_vcf/chr019");
//        countSitesinFastCallformat_fromTxt(args[0],args[1]);
            // System.out.println(Runtime.getRuntime().availableProcessors());
//        countSitesinFastCallformat_fromTxt("/Users/guoyafei/Documents/01_个人项目/05_FuGWAS/06_单倍型分析0419/test","/Users/guoyafei/Documents/01_个人项目/05_FuGWAS/06_单倍型分析0419/01_VCF/1.txt");
            //calVcfAverageDepth(args[0],args[1]);

//        BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
//        for(int i = 0; i < grid.getTaxaNumber(); i++){
//            int j;
//            for(j = 0; j < grid.getIBSDistanceMatrix().length; j++){
//                bw.write(grid.getIBSDistanceMatrix()[i][j]+" ");
//            }
//            bw.newLine();
//        }
//        bw.flush();
//        bw.close();
//        GenotypeGrid g1 = new GenotypeGrid(args[0],VCF);
//        GenotypeGrid g2 = new GenotypeGrid(args[1],VCF);
//        String infileS = args[0];
//        String outfileS = infileS.replace("vcf","ibs.txt");
//        SumTaxaDivergence std = new SumTaxaDivergence(grid);
//        std.writeDxyMatrix(outfileS, IOFileFormat.TextGzip);
//        GenotypeGrid g= GenotypeOperation.mergeGenotypesByTaxon(g1, g2);
//        SumTaxaDivergence std = new SumTaxaDivergence(g);
//        std.writeDxyMatrix(outfileS, IOFileFormat.Text);

//        grid.getIBSDistanceMatrix();
//        String[] temp=StringUtils.split("A\tB\tC");
//        for (int i = 0; i < temp.length; i++) {
//            System.out.println(temp[i]);
//        }
//        File fs = new File("/data2/yafei/003_Project3/Vmap1.1/Out/VCF/VmapE6/Alineage.all.vcf");
//        daxing.common.VCF.splitSubgenome(fs,"/data2/yafei/003_Project3/Vmap1.1/Out/VCF/VmapE6/Chr");
//        File fs1 = new File("/data2/yafei/003_Project3/Vmap1.1/Out/VCF/VmapE6/Blineage.all.vcf");
//        daxing.common.VCF.splitSubgenome(fs1,"/data2/yafei/003_Project3/Vmap1.1/Out/VCF/VmapE6/Chr");
//        File fs2 = new File("/data2/yafei/003_Project3/Vmap1.1/Out/VCF/VmapE6/Dlineage.all.vcf");
//        daxing.common.VCF.splitSubgenome(fs2,"/data2/yafei/003_Project3/Vmap1.1/Out/VCF/VmapE6/Chr");
//        daxing.common.VCF.mergeVCFtoChr("/data2/yafei/003_Project3/Vmap1.1/Out/VCF/VmapE6/Chr","/data2/yafei/003_Project3/Vmap1.1/Out/VCF/VmapE6/Merge_Chr");
//        }
//    }
}