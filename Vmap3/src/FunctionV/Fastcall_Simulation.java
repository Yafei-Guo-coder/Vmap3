package FunctionV;

import org.apache.commons.math3.distribution.GeometricDistribution;
import pgl.infra.utils.IOUtils;

import java.io.*;
import java.util.*;

public class Fastcall_Simulation {
    int haploNum = 3;
    int indiNum = 100;
    String refGenome = "/Users/guoyafei/Documents/01_个人项目/04_VmapIII/09_Fastcall2/simulation/chr1_simu.fa.gz";
    double dpi = 0.1;
    int coverage = 10;
    //30000000*10/300 reads number
    int insertL = 350;
    int readL = 150;
    String haploFile = "/Users/guoyafei/Documents/01_个人项目/04_VmapIII/09_Fastcall2/simulation/haploSet.txt";

    ArrayList<StringBuilder> haploSeqs = new ArrayList<>();
    ArrayList<ArrayList<String>> indiFa = new ArrayList<>();

    public static void main(String[] args) throws IOException {
        double mafD = 0.7;
        new Fastcall_Simulation(args);
//        GeometricDistribution distribution = new GeometricDistribution(mafD);
//        for (int i = 0; i < 100; i++) {
//            System.out.println(distribution.sample(1)[0]);
//        }
    }

    public Fastcall_Simulation(String[] args) throws IOException {
       this.haploFa(refGenome, haploNum);
//       this.simuIndi(haploSeqs, indiNum);
//       this.simuReads(indiFa, 10);
    }

    public StringBuilder writeTrueSet(int num1, int num2) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader("haploSet.txt"));
        ArrayList<String[]> haplo = null;
        String[] h = null;
        String str;
        while ((str = br.readLine()) != null){
            h = str.split("\t");
            haplo.add(h);
        }
        StringBuilder indiMut = null;

        for (int j = 0; j < haplo.size(); j++) {
            indiMut.append(haplo.get(j)[3].charAt(num1));
            indiMut.append(haplo.get(j)[3].charAt(num2));
            indiMut.append("\t");
        }
        br.close();
        return indiMut;
    }
    public StringBuilder writePosRefAlt() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader("haploSet.txt"));
        ArrayList<String[]> haplo = null;
        String[] h = null;
        String str;
        while ((str = br.readLine()) != null){
            h = str.split("\t");
            haplo.add(h);
        }
        StringBuilder indiPos = null;

        for (int m = 0; m < 3; m++) {
            for (int j = 0; j < haplo.size(); j++) {
                indiPos.append(haplo.get(j)[m]+"\t");
                indiPos.append(haplo.get(j)[m]);
            }
            indiPos.append("\n");
        }
        br.close();
        return indiPos;
    }


    public static double getRandom(double min, double max) {
        return Math.random() * (max - min) + min;
    }

    public ArrayList<StringBuilder> haploFa(String refFa, int haploNum) throws IOException {
        this.haploNum = haploNum;
        this.refGenome = refFa;
        BufferedWriter haploSet = new BufferedWriter(new FileWriter(haploFile));

        Fasta fa = new Fasta();
        fa.setBlocks(this.refGenome);
        ArrayList<ArrayList<String>> faMap = fa.getBlock();
        StringBuilder faSeq = new StringBuilder();
//        System.out.println(faMap.get(0).get(0));
//        System.out.println(faMap.get(0).get(1));
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
        BufferedWriter truebr = new BufferedWriter(new FileWriter("haploSet.txt"));
        ArrayList<StringBuilder> fas = haploFas;
        ArrayList<String> trueS = null;
        int Num = indiNum;
        ArrayList<ArrayList<String>> indiFa = null;

        for (int i = 0; i < Num; i++) {
            ArrayList<String> paired = null;
            int p1 = (int) Math.random();
            int p2 = (int) Math.random();
            paired.add(0, String.valueOf(fas.get(p1)));
            paired.add(1, String.valueOf(fas.get(p2)));
            indiFa.add(i, paired);
            trueS.add(String.valueOf(writeTrueSet(p1,p2)));
        }
        truebr.write(String.valueOf(writePosRefAlt()));
        truebr.write(String.valueOf(trueS));
        truebr.flush();
        truebr.close();
        return indiFa;
    }


    public void simuReads(ArrayList<ArrayList<String>> indiFas, int i) throws IOException {

        this.indiFa = indiFas;

        //生成reads起始位点和对应序列
        ArrayList<ArrayList<Integer>> startPos = null;
        for (int j = 0; j < indiNum; j++) {
            for (i = 0; i < coverage; i++) {
                int p = (int) getRandom(0, indiFa.get(j).get(1).length() - 351);
                startPos.get(j).add(p);
            }
        }
        ArrayList<ArrayList<String>> random350 = null;
        for (int m = 0; m < indiNum; m++) {
            for (int j = 0; j < startPos.size(); j++) {
                int ran = (int) Math.random();
                if(ran > 0.5){
                    random350.get(m).add((String) indiFa.get(0).get(1).subSequence(startPos.get(m).get(j), startPos.get(m).get(j)+350));
                }else{
                    random350.get(m).add((String) indiFa.get(0).get(2).subSequence(startPos.get(m).get(j), startPos.get(m).get(j)+350));
                }
            }
        }

        String Reads1 = null;
        String Reads2 = null;
        ArrayList<Integer> read1Q = null;
        ArrayList<Integer> read2Q = null;
        String readQ = new File("reads.Q.txt").getAbsolutePath();
        BufferedReader brQ = IOUtils.getTextReader(readQ);
        String str;
        String[] h = null;
        while ((str = brQ.readLine()) != null){
            h = str.split("\t");
            read1Q.add(Integer.valueOf(h[1]));
            read2Q.add(Integer.valueOf(h[2]));
        }
        ArrayList<Double> read1p = null;
        ArrayList<Double> read2p = null;
        for (int j = 0; j < read1Q.size(); j++) {
            double x1 = read1Q.get(j);
            double x2 = read2Q.get(j);;
            double a = 10;
            double b1 = -(x1)/10;
            double b2 = -(x2)/10;
            read1p.add(Math.pow(a, b1));
            read2p.add(Math.pow(a, b2));
        }

        for (int j = 0; j < random350.size(); j++) {
            Random r = new Random();
            String outfile1 = new File("sample"+j+"_R1.fq").getAbsolutePath();
            String outfile2 = new File("sample"+j+"_R2.fq").getAbsolutePath();
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

                bw1.write("@" + j + "\n");
                bw1.write(Reads1 + "\n");
                bw1.write("+\n");
                HashMap map = FindPos();
                int[] score = null;
                for (int l = 0; l < map.size(); l++) {
                    score[l]= (int) map.get(l);
                }

                StringBuffer sb1 = new StringBuffer();
                for (int m = 0; m < 150; m++) {
                    int in = Arrays.binarySearch(score,(int)read1Q.get(m));
                    String Q = (String) map.get(in);

//                    int score = (m - 150) * (m - 150) * 14 / 22201 + 24;
//                    String Q = getphred(score);
                    sb1.append(Q);
                }
                bw1.write(sb1.toString() + "\n");

                bw2.write("@" + j + "\n");
                bw2.write(Reads2 + "\n");
                bw2.write("+\n");
                StringBuffer sb2 = new StringBuffer();
                for (int m = 0; m < 150; m++) {
                    int in = Arrays.binarySearch(score,(int)read2Q.get(m));
                    String Q = (String) map.get(in);
//                    int score = -(m - 1) * (m - 1) * 8 / 22201 + 37;
//                    String Q = getphred(score);
                    sb2.append(Q);
                }
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