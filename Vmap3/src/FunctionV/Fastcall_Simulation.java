package FunctionV;

import org.apache.commons.math3.distribution.GeometricDistribution;
import pgl.infra.utils.IOUtils;

import java.io.*;
import java.util.*;

public class Fastcall_Simulation {
    int haploNum = 3;
    int indiNum = 2;
    String refGenome = "/Users/guoyafei/Documents/01_个人项目/04_VmapIII/09_Fastcall2/simulation/chr1_simu.fa.gz";
    String OtherGenome = "/Users/guoyafei/Documents/01_个人项目/04_VmapIII/09_Fastcall2/simulation/chr1_simu.fa.gz";
    double dpi = 0.1;
    int readsNumSingle = 4;
    //30000000*10/300 reads number

    int readsNum = (int) (readsNumSingle*0.01);
    String haploFile = "/Users/guoyafei/Documents/01_个人项目/04_VmapIII/09_Fastcall2/simulation/haploSet.txt.gz";

    ArrayList<ArrayList<String>> indiFa = new ArrayList<>();

    public static void main(String[] args) throws IOException {
        double mafD = 0.7;
        new Fastcall_Simulation(args);
    }

    public Fastcall_Simulation(String[] args) throws IOException {
        ArrayList<StringBuilder> hs = this.haploFa(refGenome, haploNum);
        ArrayList<ArrayList<String>> indi = this.simuIndi(hs, indiNum);

       this.simuReads(indi, 10);
    }



    public StringBuilder writeTrueSet(int num1, int num2) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/09_Fastcall2/simulation/haploSet.txt"));
        ArrayList<String[]> haplo = new ArrayList<>();
        String[] h = null;
        String str;
        while ((str = br.readLine()) != null){
            h = str.split("\t");
            haplo.add(h);
        }
        StringBuilder indiMut = new StringBuilder();

        for (int j = 0; j < haplo.size(); j++) {
            indiMut.append(haplo.get(j)[3].charAt(num1));
            indiMut.append(haplo.get(j)[3].charAt(num2));
            indiMut.append("\t");
        }
        br.close();
        return indiMut;
    }
    public StringBuilder writePosRefAlt() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/09_Fastcall2/simulation/haploSet.txt"));
        ArrayList<String[]> haplo = new ArrayList<>();
        String[] h = null;
        String str;
        while ((str = br.readLine()) != null){
            h = str.split("\t");
            haplo.add(h);
        }
        StringBuilder indiPos = new StringBuilder();

        for (int m = 0; m < 3; m++) {
            for (int j = 0; j < haplo.size(); j++) {
                indiPos.append(haplo.get(j)[m]+"\t");
            }
            indiPos.append("\n");
        }
        br.close();
        return indiPos;
    }


    public static double getRandom(double min, double max) {
        return Math.random() * (max - min) + min;
    }

    private ArrayList<String> simuReadsOtherChr(String otherGenome, int readsNum) throws IOException {
        this.OtherGenome = otherGenome;
        this.readsNum = readsNum;
        Fasta fa = new Fasta();
        fa.setBlocks(this.OtherGenome);
        ArrayList<ArrayList<String>> faMap = fa.getBlock();
        StringBuilder faSeq = new StringBuilder();
        faSeq.append(faMap.get(0).get(1));

        //生成reads起始位点和对应序列
        ArrayList<ArrayList<Integer>> startPos = new ArrayList<>();

        ArrayList<Integer> Pos = new ArrayList<>();

        for (int i = 0; i < readsNum; i++) {
            int p = (int) getRandom(0, faSeq.length() - 351);
            Pos.add(p);
        }

        ArrayList<ArrayList<String>> random350 = new ArrayList<>();

        ArrayList<String> se = new ArrayList<>();
        for (int j = 0; j < readsNum; j++) {
            int ran = (int) Math.random();
            if(ran > 0.5){
                int a =  Pos.get(j);
                int b = Pos.get(j)+350;
                String s = (String) faSeq.subSequence(a, b);
                se.add(j,s);
            }else{
                int a =  Pos.get(j);
                int b = Pos.get(j)+350;
                String s = (String) faSeq.subSequence(a, b);
                se.add(j,s);
            }
        }
        return se;
    }

    public ArrayList<StringBuilder> haploFa(String refFa, int haploNum) throws IOException {
        this.haploNum = haploNum;
        this.refGenome = refFa;
        BufferedWriter haploSet = new BufferedWriter(new FileWriter(haploFile));

        Fasta fa = new Fasta();
        fa.setBlocks(this.refGenome);
        ArrayList<ArrayList<String>> faMap = fa.getBlock();
        StringBuilder faSeq = new StringBuilder();
        faSeq.append(faMap.get(0).get(1));

        ArrayList<StringBuilder> haploFa = new ArrayList<>();
        for (int m = 0; m < this.haploNum; m++) {
            haploFa.add(m,faSeq);

            for (int i = 0; i < faSeq.length(); i++) {
                StringBuilder mutPos = new StringBuilder();
                char ref = faSeq.charAt(i);
                double p = getRandom(0,1);
                if (p <= dpi) {
                    char alt = randomAllele(faSeq.charAt(i));
                    mutPos.append(i+"\t"+ref+"\t"+alt+"\t");
                    Double maf = null;
                    double mafD = 0.7;
                    GeometricDistribution distribution = new GeometricDistribution(mafD);
                    int D = distribution.sample(1)[0];
                    double maf1 = 0;
                    double p2 = 0;
                    if(D==0){
                        maf1 = getRandom(0,0.1);
                        p2 = getRandom(0,1);
                    }else if(D==1){
                        maf1 = getRandom(0.1,0.2);
                        p2 = getRandom(0,1);
                    }else if(D==2){
                        maf1 = getRandom(0.2,0.3);
                        p2 = getRandom(0,1);
                    }else if(D==3){
                        maf1 = getRandom(0.3,0.4);
                        p2 = getRandom(0,1);
                    }else if(D==4){
                        maf1 = getRandom(0.4,0.5);
                        p2 = getRandom(0,1);
                    }else{
                        maf1 = getRandom(0,0.1);
                        p2 = getRandom(0,1);
                    }
                    for (int q = 0; q < haploNum; q++) {
                        if (maf1 > p2) {
                            haploFa.get(m).setCharAt(i, alt);
                            mutPos.append(alt);
                        }else{
                            mutPos.append(ref);
                        }
                    }
                    haploSet.write(mutPos.toString());
                    haploSet.newLine();
                }
            }
        }
        haploSet.flush();
        haploSet.close();
        return haploFa;
    }

    public ArrayList<ArrayList<String>> simuIndi(ArrayList<StringBuilder> haploFas, int indiNum) throws IOException {
        BufferedWriter truebr = new BufferedWriter(new FileWriter("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/09_Fastcall2/simulation/TrueSet1.txt"));
        ArrayList<StringBuilder> fas = haploFas;
//        ArrayList<String> trueS = new ArrayList<>();
        int Num = indiNum;
        ArrayList<ArrayList<String>> indiFa = new ArrayList<>();
        truebr.write(String.valueOf(writePosRefAlt()));
        for (int i = 0; i < Num; i++) {
            ArrayList<String> paired = new ArrayList<>();
            int p1 = (int) Math.random();
            int p2 = (int) Math.random();
            paired.add(0, String.valueOf(fas.get(p1)));
            paired.add(1, String.valueOf(fas.get(p2)));
            indiFa.add(i, paired);
            truebr.write(String.valueOf(writeTrueSet(p1,p2)));
            truebr.write("\n");
        }
//        truebr.write(String.valueOf(trueS));

        truebr.flush();
        truebr.close();
        return indiFa;
    }
    public void simuReads(ArrayList<ArrayList<String>> indiFas, int i) throws IOException {
        ArrayList<String> se2 = this.simuReadsOtherChr(OtherGenome, readsNum);

        this.indiFa = indiFas;

        //生成reads起始位点和对应序列
        ArrayList<ArrayList<Integer>> startPos = new ArrayList<>();
        for (int j = 0; j < indiNum; j++) {
            ArrayList<Integer> Pos = new ArrayList<>();
            for (i = 0; i < readsNumSingle; i++) {
                int p = (int) getRandom(0, indiFa.get(j).get(1).length() - 351);
                Pos.add(p);
            }
            startPos.add(j,Pos);
        }
        ArrayList<ArrayList<String>> random350 = new ArrayList<>();
        for (int m = 0; m < indiNum; m++) {
            ArrayList<String> se = new ArrayList<>();
            for (int j = 0; j < readsNumSingle; j++) {
                int ran = (int) Math.random();
                if(ran > 0.5){
                    int a =  startPos.get(m).get(j);
                    int b = startPos.get(m).get(j)+350;
                    String s = (String) indiFa.get(m).get(1).subSequence(a, b);
                    se.add(j,s);
                }else{
                    int a =  startPos.get(m).get(j);
                    int b = startPos.get(m).get(j)+350;
                    String s = (String) indiFa.get(m).get(1).subSequence(a, b);
                    se.add(j,s);
                }
            }
            random350.add(m,se);
        }
        String Reads1 = new String();
        String Reads2 = new String();
        ArrayList<Integer> read1Q = new ArrayList<>();
        ArrayList<Integer> read2Q = new ArrayList<>();
        String readQ = new File("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/09_Fastcall2/simulation/reads.Q.txt").getAbsolutePath();
        BufferedReader brQ = IOUtils.getTextReader(readQ);
        String str;
        String[] h = null;
        int count = 0;
        while ((str = brQ.readLine()) != null){
            h = str.split("\t");
            read1Q.add(count, (int)Double.parseDouble(h[1]));
            read2Q.add(count, (int)Double.parseDouble(h[2]));
            count++;
        }
        ArrayList<Double> read1p = new ArrayList<>();
        ArrayList<Double> read2p = new ArrayList<>();
        for (int j = 0; j < read1Q.size(); j++) {
            double x1 = read1Q.get(j);
            double x2 = read2Q.get(j);;
            double a = 10;
            double b1 = -(x1)/10;
            double b2 = -(x2)/10;
            read1p.add(Math.pow(a, b1));
            read2p.add(Math.pow(a, b2));
        }

        HashMap map = FindPos();
        StringBuffer sb2 = new StringBuffer();
        for (int m = 0; m < 150; m++) {
            String Q = (String) map.get(read2Q.get(m));
            sb2.append(Q);
        }
        StringBuffer sb1 = new StringBuffer();
        for (int m = 0; m < 150; m++) {
            String Q = (String) map.get(read1Q.get(m));
            sb1.append(Q);
        }
        for (int j = 0; j < random350.size(); j++) {
            Random r = new Random();
            String outfile1 = new File("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/09_Fastcall2/simulation/sample"+j+"_R1.fq").getAbsolutePath();
            String outfile2 = new File("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/09_Fastcall2/simulation/sample"+j+"_R2.fq").getAbsolutePath();
            BufferedWriter bw1 = IOUtils.getTextWriter(outfile1);
            BufferedWriter bw2 = IOUtils.getTextWriter(outfile2);
            for (int k = 0; k < random350.get(j).size(); k++) {
                Reads1 = random350.get(j).get(k).substring(0, 150);
                StringBuilder sbreads1 = new StringBuilder();
                for (int m = 0; m < Reads1.length(); m++) {

                    String index = Reads1.substring(m, m + 1);
                    String index1 = null;
                    double ran = getRandom(0,1);
                    if(ran < read1p.get(m)){
                        index1 = String.valueOf(randomAllele(Reads1.charAt(m)));
                        sbreads1.append(index1);
                    } else {
                        sbreads1.append(index);
                    }
                }
                Reads1 = sbreads1.toString();
                Reads2 = random350.get(j).get(k).substring(200, 350);
                Reads2 = new StringBuffer(Reads2).reverse().toString();
                StringBuilder sbreads2 = new StringBuilder();
                for (int l = 0; l < Reads2.length(); l++) {
                    String index = Reads2.substring(l, l + 1);
                    String index1 = null;
                    double ran = getRandom(0,1);
                    if (ran < read2p.get(l)) {
                        index1 = String.valueOf(randomAllele(Reads2.charAt(l)));
                        sbreads2.append(index1);
                    } else {
                        sbreads2.append(index);
                    }
                }
                Reads2 = sbreads2.toString();

                StringBuffer sb = new StringBuffer();
                for (int m = 0; m < Reads2.length(); m++) {
                    String index = Reads2.substring(m, m + 1);
                    String index1 = null;
                    if (index.equals("A")) {
                        index1 = "T";
                    }
                    if (index.equals("G")) {
                        index1 = "C";
                    }
                    if (index.equals("T")) {
                        index1 = "A";
                    }
                    if (index.equals("C")) {
                        index1 = "G";
                    }
                    sb.append(index1);
                }
                Reads2 = sb.toString();

                bw1.write("@" + k + "\n");
                bw1.write(Reads1 + "\n");
                bw1.write("+\n");
                bw1.write(sb1.toString() + "\n");

                bw2.write("@" + k + "\n");
                bw2.write(Reads2 + "\n");
                bw2.write("+\n");
                bw2.write(sb2.toString() + "\n");
            }
            for (int k = 0; k < se2.size(); k++) {
                Reads1 = se2.get(k).substring(0, 150);
                StringBuilder sbreads1 = new StringBuilder();
                for (int m = 0; m < Reads1.length(); m++) {

                    String index = Reads1.substring(m, m + 1);
                    String index1 = null;
                    double ran = getRandom(0,1);
                    if(ran < read1p.get(m)){
                        index1 = String.valueOf(randomAllele(Reads1.charAt(m)));
                        sbreads1.append(index1);
                    } else {
                        sbreads1.append(index);
                    }
                }
                Reads1 = sbreads1.toString();
                Reads2 = se2.get(k).substring(200, 350);
                Reads2 = new StringBuffer(Reads2).reverse().toString();
                StringBuilder sbreads2 = new StringBuilder();
                for (int l = 0; l < Reads2.length(); l++) {
                    String index = Reads2.substring(l, l + 1);
                    String index1 = null;
                    double ran = getRandom(0,1);
                    if (ran < read2p.get(l)) {
                        index1 = String.valueOf(randomAllele(Reads2.charAt(l)));
                        sbreads2.append(index1);
                    } else {
                        sbreads2.append(index);
                    }
                }
                Reads2 = sbreads2.toString();

                StringBuffer sb = new StringBuffer();
                for (int m = 0; m < Reads2.length(); m++) {
                    String index = Reads2.substring(m, m + 1);
                    String index1 = null;
                    if (index.equals("A")) {
                        index1 = "T";
                    }
                    if (index.equals("G")) {
                        index1 = "C";
                    }
                    if (index.equals("T")) {
                        index1 = "A";
                    }
                    if (index.equals("C")) {
                        index1 = "G";
                    }
                    sb.append(index1);
                }
                Reads2 = sb.toString();

                bw1.write("@" + k + "_other\n");
                bw1.write(Reads1 + "\n");
                bw1.write("+\n");
                bw1.write(sb1.toString() + "\n");

                bw2.write("@" + k + "_other\n");
                bw2.write(Reads2 + "\n");
                bw2.write("+\n");
                bw2.write(sb2.toString() + "\n");
            }

            bw1.flush();
            bw1.close();
            bw2.flush();
            bw2.close();
            }
        }

    public static char randomAllele(char alle){
        char alt = 0;
        Character a = alle;
        double p = getRandom(0,1);

        HashMap map = new HashMap<Integer, Character>();
        map.put(0, 'A');
        map.put(1, 'T');
        map.put(2, 'C');
        map.put(3, 'G');

        int in = 0;
        for (int j = 0; j < 4; j++) {
            if(a.equals(map.get(j))){
                in = j;
            }
        }
        StringBuilder sb = new StringBuilder();
        for (int m = 0; m < 4; m++) {
            if(m==in)continue;
            sb.append(map.get(m));
        }

        if(p <= 0.33){
            alt = sb.charAt(0);
        }else if(0.33 < p && p <= 0.66){
            alt = sb.charAt(1);
        }else{
            alt = sb.charAt(2);
        }

        return alt;
    }
    public HashMap FindPos(){
        HashMap map = new HashMap<Integer, String>();
        List temp = new ArrayList<String>();
        map.put(0, "!");
        map.put(1, "\"");
        map.put(2, "#");
        map.put(3, "$");
        map.put(4, "%");
        map.put(5, "&");
        map.put(6, "'");
        map.put(7, "(");
        map.put(8, ")");
        map.put(9, "*");
        map.put(10, "+");
        map.put(11, ",");
        map.put(12, "-");
        map.put(13, ".");
        map.put(14, "/");
        map.put(15, "0");
        map.put(16, "1");
        map.put(17, "2");
        map.put(18, "3");
        map.put(19, "4");
        map.put(20, "5");
        map.put(21, "6");
        map.put(22, "7");
        map.put(23, "8");
        map.put(24, "9");
        map.put(25, ":");
        map.put(26, ";");
        map.put(27, "<");
        map.put(28, "=");
        map.put(29, ">");
        map.put(30, "?");
        map.put(31, "@");
        map.put(32, "A");
        map.put(33, "B");
        map.put(34, "C");
        map.put(35, "D");
        map.put(36, "E");
        map.put(37, "F");
        map.put(38, "G");
        map.put(39, "H");
        map.put(40, "I");
        return map;
    }

}