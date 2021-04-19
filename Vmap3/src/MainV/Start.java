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
//import static FunctionV.CountSite.countrepeatIndelinFastCallformat;

public class Start {
    public static void main(String[] args) throws IOException {
        countSitesinFastCallformat_fromTxt(args[0]);
    }
}
