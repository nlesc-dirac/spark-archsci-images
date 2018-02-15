public class Driver
{
  public native String stringMethod(String text);
  public native int readMSAndBack(String[] argc);


  public static void main(String[] args)
  {
    System.loadLibrary("Driver");
    //System.load("Driver");
    Driver driver=new Driver();

   // String text=driver.stringMethod("sm1.ms");
   // System.out.println("stringMethod: "+text);
   
    int I1=driver.readMSAndBack(args);
    System.out.println("quitting with output "+I1);
  }
}
