/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.data4CandChIA_PET;

import java.io.BufferedReader;
import java.io.BufferedWriter;

/**
 *
 * @author xuebozhao
 */
public class GetCircosBed {
    
    GetCircosBed(String infileS,String outfileS) {
        
        this.sum(infileS,outfileS);
        
    }
    
    public void readBed (String infileS,String outfileS){
        try{
            
            BufferedReader br;            
            br = XueboIOUtils.getTextReader(infileS);
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            
            String temp = null;
            String L1 = null;
            String L2 = null;
            int i = 0;

            while (( temp = br.readLine()) != null) {
                
                    ++i;
                    
                    if (i % 100 == 0) {

                        System.out.println("Filtering " + i + "....");

                    }

                    String[] tem = temp.split("\t");
                    
                    int aa = Integer.valueOf(tem[2]) + 10 ;
                    
                    L1 = "chr" + tem[0] + "\t" + tem[1] + "\t" + aa + "\t" + "chr2" + "\t" + "7808877" + "\t" + "7810811";
                    
                    bw.write(L1 + "\n"); 
            }
            
            bw.flush(); 
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
    
    
    public void sum (String infileS,String outfileS){
        try{
            
            BufferedReader br;            
            br = XueboIOUtils.getTextReader(infileS);
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            
            String temp = null;
            int L1 = 0;
            int L2 = 0;
            int i = 0;

            while (( temp = br.readLine()) != null) {
                
                    ++i;
                    
                    if (i % 2 == 0) {

                        System.out.println("Filtering " + i + "....");

                    }

                    if (temp.startsWith("CHROM"))
                        continue;
                    String[] tem = temp.split("\t");
                    
                    int aa = Integer.valueOf(tem[1]);
                    
                    L1 = L1 + aa;
                    
                    
                    //bw.write(L1 + "\n"); 
            }
            bw.write(L1 + "\n");
            bw.flush(); 
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
