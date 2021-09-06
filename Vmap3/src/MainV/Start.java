package MainV;


import it.unimi.dsi.fastutil.ints.IntArrayList;
import pgl.app.fastCall2.IndividualGenotype;
import pgl.infra.utils.IOUtils;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class Start {
    public static void main(String[] args) throws IOException {
        String[] taxaNames = null;
        String inputDirS = "/Users/guoyafei/Documents/01_个人项目/04_VmapIII/09_Fastcall2/ing";
        String outputDirS = "/Users/guoyafei/Documents/01_个人项目/04_VmapIII/09_Fastcall2/ing2";
        File outDir = new File (outputDirS);
//        outDir.mkdir();

        List<File> chrDirList = IOUtils.getDirListInDir(inputDirS);
        File[] chrOutDirs = new File[chrDirList.size()];
        for (int i = 0; i < chrDirList.size(); i++) {
            File f = new File (outDir, chrDirList.get(i).getName());
            f.mkdir();
            chrOutDirs[i] = f;
        }
        for (int i = 0; i < chrDirList.size(); i++) {
            List<File> taxaDirlist = IOUtils.getDirListInDir(String.valueOf(chrDirList.get(i)));
            taxaNames = new String[taxaDirlist.size()];
            for (int j = 0; j < taxaDirlist.size(); j++) {
                File f = new File(chrOutDirs[i], taxaDirlist.get(j).getName());
                f.mkdir();
            }
        }
        for (int i = 0; i < chrDirList.size(); i++) {
            List<File> taxaDirlist = IOUtils.getDirListInDir(String.valueOf(chrDirList.get(i)));
            taxaNames = new String[taxaDirlist.size()];
            for (int j = 0; j < taxaDirlist.size(); j++) {
                taxaNames[j] = taxaDirlist.get(j).getName();
            }
            Arrays.sort(taxaNames);
            for (int k = 0; k < taxaDirlist.size(); k++) {
//                String taxaS = new File(String.valueOf(taxaDirlist.get(k))).getAbsolutePath();
//                String outtaxaS = taxaS.replace("ing", "ing2");
//                File f = new File (outtaxaS);
//                f.mkdir();
                List<File> ingList = Arrays.asList(taxaDirlist.get(k).listFiles());
                for (int j = 0; j < ingList.size(); j++) {
                    String fileS = new File(String.valueOf(ingList.get(j))).getAbsolutePath();
//                        System.out.println(fileS);
//                    IndividualGenotype ing = new IndividualGenotype(fileS);
                    DataInputStream dis = IOUtils.getBinaryGzipReader(fileS);
                    String taxonName = dis.readUTF();
                    short chrom = dis.readShort();
                    int binStart = dis.readInt();
                    int binEnd = dis.readInt();
                    IntArrayList codedAlleleInfo = null;
                    codedAlleleInfo = new IntArrayList();
                    int currentRecord = 0;
                    while ((currentRecord = dis.readInt()) != Integer.MIN_VALUE) {
                        codedAlleleInfo.add(currentRecord);
                    }
                    for (int m = 0; m < codedAlleleInfo.size(); m++) {
                        int v = codedAlleleInfo.getInt(m) & 192;
                        if (v != 192) continue;
//                            System.out.println(Integer.toBinaryString(codedAlleleInfo.getInt(m)));
//                            System.out.println(ing.getAlleleBase(m));
//                            v = codedAlleleInfo.getInt(m) & 2147483392;
//                            v = v + 191;
                        codedAlleleInfo.set(m,codedAlleleInfo.getInt(m)-1);
//                            ing.codedAlleleInfo.set(m,codedAlleleInfo.getInt(m)-1);

//                            System.out.println(Integer.toBinaryString(codedAlleleInfo.getInt(m)));
//                            System.out.println(ing.getAlleleBase(m));
//                            System.out.println("------------------");
                    }
                    String sb = new File(String.valueOf(ingList.get(j))).getAbsolutePath();
                    String sb2 = sb.replaceFirst("ing", "ing2");
                    DataOutputStream dos = IOUtils.getBinaryGzipWriter(sb2);
                    try {
                        dos.writeUTF(taxonName);
                        dos.writeShort((short)chrom);
                        dos.writeInt(binStart);
                        dos.writeInt(binEnd);
                    }
                    catch (Exception e) {
                        e.printStackTrace();
                        System.exit(1);
                    }
                    for (int m = 0; m < codedAlleleInfo.size(); m++) {
                        dos.writeInt(codedAlleleInfo.getInt(m));
                    }
                    dos.flush();
                    dos.close();
                }
            }
//        splitVcf(args[0]);
            //splitVcf("/data2/yafei/004_Vmap3/VCF/Raw_VCF/AA_vcf/chr019");
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
        }
    }
}