package MainV;

import org.apache.commons.cli.*;
import pgl.AppNames;
import pgl.PGLAPPEntrance;
import pgl.app.fastCall.FastCall;
import pgl.app.fastCall2.FastCall2;
import pgl.app.hapScanner.HapScanner;
import pgl.app.popdep.PopDep;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.summa.SumTaxaDivergence;
import pgl.infra.utils.CLIInterface;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

import static FunctionV.CountSite.*;
import static pgl.infra.dna.genot.GenoIOFormat.VCF;
import static pgl.infra.dna.genot.GenoIOFormat.VCF_GZ;
//import static FunctionV.CountSite.countrepeatIndelinFastCallformat;

public class Start {
    public static void main(String[] args) throws IOException {
//        splitVcf(args[0]);
        //splitVcf("/data2/yafei/004_Vmap3/VCF/Raw_VCF/AA_vcf/chr019");
//        countSitesinFastCallformat_fromTxt(args[0],args[1]);
       // System.out.println(Runtime.getRuntime().availableProcessors());
//        countSitesinFastCallformat_fromTxt("/Users/guoyafei/Documents/01_个人项目/05_FuGWAS/06_单倍型分析0419/test","/Users/guoyafei/Documents/01_个人项目/05_FuGWAS/06_单倍型分析0419/01_VCF/1.txt");
        calVcfAverageDepth(args[0],args[1]);
//        GenotypeGrid grid = new GenotypeGrid(args[0],VCF_GZ);
//        String infileS = args[0];
//        String outfileS = infileS.replace("vcf.gz","ibs.txt.gz");
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
//        SumTaxaDivergence std = new SumTaxaDivergence(grid);
//        std.writeDxyMatrix(outfileS, IOFileFormat.TextGzip);
        //grid.getIBSDistanceMatrix();
    }
}
