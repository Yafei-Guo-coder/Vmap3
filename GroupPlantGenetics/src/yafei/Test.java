package yafei;

import gnu.trove.list.array.TIntArrayList;

public class Test {

    public static void main(String[] args) {
        TIntArrayList tIntArrayList=  new TIntArrayList();
        tIntArrayList.add(5);
        tIntArrayList.add(6);
        tIntArrayList.add(7);
        for (int i = 0; i < tIntArrayList.size(); i++) {
            System.out.println(tIntArrayList.get(i));
        }
        System.out.println("ok");
    }
}
