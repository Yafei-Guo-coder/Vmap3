package MainV;

import org.apache.commons.cli.*;
import pgl.AppNames;
import pgl.PGLAPPEntrance;
import pgl.app.fastCall.FastCall;
import pgl.app.fastCall2.FastCall2;
import pgl.app.hapScanner.HapScanner;
import pgl.app.popdep.PopDep;
import pgl.infra.utils.CLIInterface;

import java.io.File;
import java.io.IOException;

import static FunctionV.CountSite.countSitesinFastCallformat_fromTxt;
import static FunctionV.CountSite.splitVcf;
//import static FunctionV.CountSite.countrepeatIndelinFastCallformat;

public class Start {
    public static void main(String[] args) throws IOException {
        splitVcf(args[0]);
       // splitVcf("/Users/guoyafei/Desktop/Test");
//        countSitesinFastCallformat_fromTxt(args[0],args[1]);
       // System.out.println(Runtime.getRuntime().availableProcessors());
//        countSitesinFastCallformat_fromTxt("/Users/guoyafei/Documents/01_个人项目/05_FuGWAS/06_单倍型分析0419/test","/Users/guoyafei/Documents/01_个人项目/05_FuGWAS/06_单倍型分析0419/01_VCF/1.txt");
    }
}
