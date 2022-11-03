import java.util.*;
public class Damtpart1 {
    static void sampleGen(int numberofsamples,String distribution,double p)
    {
        //checks condition and calls bernouli and geometric functions based on user input
        if(distribution.toLowerCase().equals("bernoulli")){
            bernoulisamplegen(numberofsamples,distribution,p);
        }
        if(distribution.toLowerCase().equals("geometric")){
            geometricsamplegen(numberofsamples,distribution,p);
        }
    }
    static void sampleGen(int numberofsamples,String distribution,int n,double p)
    {
        //checks condition and calls binomial and negbinomial functions based on user input
        if(distribution.toLowerCase().equals("binomial")){
            binomialsamplegen(numberofsamples,distribution,n,p);
        }
        if(distribution.toLowerCase().equals("neg-binomial")){
            negbinomialsamplegen(numberofsamples,distribution,n,p);
        }
    }
    static void sampleGen(int numberofsamples,String distribution,int alpha,int lambda)
    {
        //checks condition and calls gamma and uniform functions based on user input
        if(distribution.toLowerCase().equals("gamma")){
            gammasamplegen(numberofsamples,distribution,alpha,lambda);
        }
        if(distribution.toLowerCase().equals("uniform")){
            uniformsamplegen(numberofsamples,distribution,alpha,lambda);
        }
    }
    static void sampleGen(int numberofsamples,String distribution,int lambda)
    {
        //checks conditions and calls exponential and poisson functions based on user input
        if(distribution.toLowerCase().equals("exponential")){
            exponentialsamplegen(numberofsamples,distribution,lambda);
        }
        if(distribution.toLowerCase().equals("poisson")){
            poissonsamplegen(numberofsamples,distribution,lambda);
        }
    }

    static void sampleGen(int numberofsamples,String distribution,double a,double b){
        if(distribution.toLowerCase().equals("normal")){
            normalsamplegen(numberofsamples,distribution,a,b);
        }
    }
    //Random number generator for bernouli distributions
    static void bernoulisamplegen(int numberofsamples,String distribution,double p)
    {
        //Seed initialized to 5
        Random r=new Random(5);
        //used to store samples
        int output[]=new int[numberofsamples];
        int index=0;
        for(int i=0;i<numberofsamples;i++){
            double d=r.nextDouble();
            if(d<p){
                output[index]=1;
                index++;
            }
            else{
                output[index]=0;
                index++;
            }
        }
        //Prints the output array
        printArray(output);
    }

    //Random number generator for binomial distribution
    static void binomialsamplegen(int numberofsamples,String distribution,int n,double p)
    {
        //Seed initialized to 5
        Random r=new Random(5);
        //Output array that stores the samples
        int output[]=new int[numberofsamples];
        int count=0;
        int index=0;
		for(int i=0;i<numberofsamples;i++){
		    count=0;
		    for(int j=0;j<n;j++){
		        double d=r.nextDouble();
		        if(d<p){
		            count++;
		        }
		    }
		    output[index]=count;
            index++;
		}
        //Prints the output array
        printArray(output);
    }

    //Random number generator for geometric distribution
    static void geometricsamplegen(int numberofsamples,String distribution,double p)
    {
        Random r=new Random(5);
        int count=0;
        //Output array that stores the samples
        int output[]=new int[numberofsamples];
        int index=0;
        double d=r.nextDouble();
        for(int i=0;i<numberofsamples;i++)
        {
            if(i>0)
            {
                d=r.nextDouble();
            }
            count=0;
            while(d>p)
            {
                count++;
                d=r.nextDouble();
            }
            count=count+1;
            output[index]=count;
            index++;
        }
        //Prints the output array
        printArray(output);
    }

    //Random number generator for Negative Binomial distribution
    static void negbinomialsamplegen(int numberofsamples,String distribution,int k,double p)
    {
        //Seed initialized to 5
        Random r=new Random(5);
        int count=0;
		int num=0;
		double d;
        //Output array to store the samples
        int output[]=new int[numberofsamples];
        int index=0;
	    for(int i=0;i<numberofsamples;i++){
            count=0;
            num=0;
            if(i>0)
            {
                d=r.nextDouble();
                if(d<p)
                {
                    count++;
                    num++;
                }
                else
                {
                    count++;
                }
            }
            //runs until it find k-successes
            while(num!=k)
            {
                d=r.nextDouble();
                if(d<p){
                    count++;
                    num++;
                }else{
                    count++;
                }
	        }
            output[index]=count;
            index++;
		}
        //Print the output array that generates samples
        printArray(output);
    }

    //Random Number generator for poisson distribution
    static void poissonsamplegen(int numberofsamples,String distribution,int lambda)
    {
        //Seed initialized to 5
        Random r=new Random(5);
        double d;
        int count=0;
        double func=0;
        //Output array that stores the samples
        int output[]=new int[numberofsamples];
        int index=0;
        for(int i=0;i<numberofsamples;i++){
            d=r.nextDouble();
            func=Math.exp(-lambda);
            count=0;
            while(d>=func){
                func=func+Math.exp(-lambda)*Math.pow(lambda,count)/factorial(count);
                count++;
            }
            count=count-1;
            output[index]=count;
            index++;
        }
        //Printing the output array
        printArray(output);
    }

    //Random number generator for uniform distribution
    static void uniformsamplegen(int numberofsamples,String dist,int a,int b)
    {
        //Seed initialized to 5
        Random r=new Random(5);
        double d;
        double x=0;
        double output[]=new double[numberofsamples];
        int index=0;
        for(int i=0;i<numberofsamples;i++){
            d=r.nextDouble();
            x=a+d*(b-a);
            output[index]=x;
            index++;
        }
        //Prints the output array
        printArray(output);
    }

    //Random number generator for exponential distribution
    static void exponentialsamplegen(int numberofsamples,String distribution,int lambda)
    {
        Random r=new Random(5);
        double l=lambda;
        double output[]=new double[numberofsamples];
        int index=0;
        for(int i=0;i<numberofsamples;i++){
            double d=r.nextDouble();
            double sample=1/l;
            double mathcal=Math.log(1-d);
            double cal=-(sample*mathcal);
            output[index]=cal;
            index++;
        }
        //Prints the output array
        printArray(output);
    }

    //Random number generator for Gamma distribution
    static void gammasamplegen(int numberofsamples,String distribution,int alpha,int lambda)
    {
        //Seed initialized to 5
        Random r=new Random(5);
        double l=lambda;
        double sum=0;
        double sample,mathcal,cal;
        double output[]=new double[numberofsamples];
        int index=0;
        for(int i=0;i<numberofsamples;i++)
        {
            for(int j=0;j<alpha;j++)
            {
                double d=r.nextDouble();
                sample=1/l;
                mathcal=Math.log(1-d);
                cal=-(sample*mathcal);
                sum=sum+cal;
            }
            output[index]=sum;
            index++;
        }
        //Prints the output array
        printArray(output);
    }

    static void normalsamplegen(int numofsamples,String distribution,double a,double b){
        double samples[] = new double[numofsamples];
        int i=0;
        //Seed Initialization 1
        Random temp1=new Random(5);
        //Seed Initialization 2
        Random temp2=new Random(5);
        while(i<(numofsamples/2)+1){
            double s=temp1.nextDouble();
            double q=temp2.nextDouble();
            float c = (float) Math.sqrt(-2*Math.log(s));
            float d = (float) (2*Math.PI*q);
            double val1 = (double)c*(Math.cos(d));
            double val2 = c*Math.sin(d);
            samples[i++] = val1*b + a;
            samples[i++] = val2*b+a;
        }
        //Function to print sample array
        printArray(samples);
    }
    //Function to print integer array
    static void printArray(int outputvalues[])
    {
        for(int i=0;i<outputvalues.length;i++){
            System.out.print(outputvalues[i]+" ");
        }
    }

    //Function to print array of double type
    static void printArray(double outputvalues[])
    {
        for(int i=0;i<outputvalues.length;i++){
            System.out.print(outputvalues[i]+" ");
        }
    }
    
    //Function that return factorial of a number-used in poisson distribution
    static int factorial(int n)
    {
        int f=1;
        for(int i=1;i<=n;i++){
            f=f*i;
        }
        return f;
    }

    //Main Method starts here
    public static void main(String args[])
    {
        //Program starts here
        //Command Line arguments to take samples,distribution and parameters
        int numberofsamples=Integer.parseInt(args[0]);
        String distribution=args[1];
        //conditions to check type of distribution and invoke parameters based on distribution given
        if(distribution.toLowerCase().equals("binomial")||distribution.toLowerCase().equals("neg-binomial"))
        {
            int n=Integer.parseInt(args[2]);
            double p=Float.parseFloat(args[3]);
            sampleGen(numberofsamples,distribution,n,p);
        }
        else if(distribution.toLowerCase().equals("bernoulli"))
        {
            double p=Float.parseFloat(args[2]);
            sampleGen(numberofsamples,distribution,p);
        }
        else if(distribution.toLowerCase().equals("geometric"))
        {
            double lambda=Double.parseDouble(args[2]);
            sampleGen(numberofsamples,distribution,lambda);
        }
        else if(distribution.toLowerCase().equals("exponential")||distribution.toLowerCase().equals("poisson"))
        {
            int lambda=Integer.parseInt(args[2]);
            sampleGen(numberofsamples,distribution,lambda);
        }
        else if(distribution.toLowerCase().equals("gamma")||distribution.toLowerCase().equals("uniform"))
        {
            int alpha=Integer.parseInt(args[2]);
            int lambda=Integer.parseInt(args[3]);
            sampleGen(numberofsamples,distribution,alpha,lambda);
        }
        else if(distribution.toLowerCase().equals("normal")){
            double mean=Double.parseDouble(args[2]);
            double std=Double.parseDouble(args[3]);
            sampleGen(numberofsamples,distribution,mean,std);
        }
    }
 
}
