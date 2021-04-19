package MainV;

import java.io.IOException;

import static FunctionV.CountSite.countSitesinFastCallformat_fromTxt;
//import static FunctionV.CountSite.countrepeatIndelinFastCallformat;

public class Start {
    public static void main(String[] args) throws IOException {
        //countrepeatIndelinFastCallformat("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Test_chr001/05_VcfCheck");
        countSitesinFastCallformat_fromTxt("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Test_chr001/05_VcfCheck");
    }
}