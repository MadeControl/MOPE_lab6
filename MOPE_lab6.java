import java.util.*;

public class MOPE_lab6 {

    public static int minX1 = 20;
    public static int maxX1 = 70;
    public static int minX2 = -15;
    public static int maxX2 = 45;
    public static int minX3 = 20;
    public static int maxX3 = 35;
    public static int middleMinX;
    public static int middleMaxX;
    public static double minY;
    public static double maxY;
    public static double x01;
    public static double x02;
    public static double x03;
    public static double deltaX1;
    public static double deltaX2;
    public static double deltaX3;
    public static double[] averageY;
    private static Random r = new Random();
    private static double[][] matrixX;
    private static Double[][] normMatrix;
    private static int N = 15;
    private static GaussianElimination solveMatrix = new GaussianElimination();
    private static Data data = new Data();
    private static List<Double> dispersionY = new ArrayList<Double>();
    private static boolean bool = false;
    private static double dispersionB2;
    private static int f3;
    private static double[] stY;
    private static int coef = 0;
    private static int nL = 0;
    private static int nEf = 0;
    private static int nQt = 0;

    private static int m = 2;
    public static int d = 0;
    private static double p = 0.95;
    private static double q = 0.05;

    public static void main(String[] args) {
        for(int i=0; i < 100; i++) {
            count();
            System.out.println();
        }
        System.out.println("Кількість рівнянь лінійної форми: "+nL+"\nКількість рівнянь з ефектом взаємодії: "
                +nEf+"\nКількість рівнянь з квадратичними членами: "+nQt);
    }

    public static boolean count(){

        System.out.println("------------------------------------------------------------------------------------------");

        middleMinX = (minX1 + minX2 + minX3)/3;
        middleMaxX = (maxX1 + maxX2 + maxX3)/3;

        minY = 200 + middleMinX;
        maxY = 200 + middleMaxX;

        x01 = (minX1 + maxX1) / 2.;
        x02 = (minX2 + maxX2) / 2.;
        x03 = (minX3 + maxX3) / 2.;

        deltaX1 = maxX1 - x01;
        deltaX2 = maxX2 - x02;
        deltaX3 = maxX3 - x03;

        normMatrix = new Double[][]{
                {-1., -1., -1., +1., +1., +1., -1., +1., +1., +1.},
                {-1., -1., +1., +1., -1., -1., +1., +1., +1., +1.},
                {-1., +1., -1., -1., +1., -1., +1., +1., +1., +1.},
                {-1., +1., +1., -1., -1., +1., -1., +1., +1., +1.},
                {+1., -1., -1., -1., -1., +1., +1., +1., +1., +1.},
                {+1., -1., +1., -1., +1., -1., -1., +1., +1., +1.},
                {+1., +1., -1., +1., -1., -1., -1., +1., +1., +1.},
                {+1., +1., +1., +1., +1., +1., +1., +1., +1., +1.},
                {-1.215, 0., 0., 0., 0., 0., 0., 1.4623, 0., 0.},
                {+1.215, 0., 0., 0., 0., 0., 0., 1.4623, 0., 0.},
                {0., -1.215, 0., 0., 0., 0., 0., 0., 1.4623, 0.},
                {0., +1.215, 0., 0., 0., 0., 0., 0., 1.4623, 0.},
                {0., 0., -1.215, 0., 0., 0., 0., 0., 0., 1.4623},
                {0., 0., +1.215, 0., 0., 0., 0., 0., 0., 1.4623},
                {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}
        };

        matrixX = new double[15][10];

        double x1=0,x2=0,x3=0, xLst[];
        for(int i = 0; i < matrixX.length; i++){
            if(i < 8){
                x1 = normMatrix[i][0] == 1. ? maxX1 : minX1;
                x2 = normMatrix[i][1] == 1. ? maxX2 : minX2;
                x3 = normMatrix[i][2] == 1. ? maxX3 : minX3;
            }
            else{
                xLst = countX(normMatrix[i][0], normMatrix[i][1], normMatrix[i][2]);
                x1 = xLst[0];
                x2 = xLst[1];
                x3 = xLst[2];
            }
            matrixX[i] = new double[]{x1, x2, x3, x1*x2, x1*x3, x2*x3, x1*x2*x3, Math.pow(x1,2), Math.pow(x2,2), Math.pow(x3,2)};
        }

        double[][] matrixY = generateMatrix(N,m);

        averageY = getAverage(matrixY,1);

        double[] mx_i = getAverage(matrixX, 0);
        double my = sum(averageY)/15;

        double[][] unknown = new double[][]{
                {1., mx_i[0], mx_i[1], mx_i[2], mx_i[3], mx_i[4], mx_i[5], mx_i[6], mx_i[7], mx_i[8], mx_i[9]},
                {mx_i[0], a(1, 1), a(1, 2), a(1, 3), a(1, 4), a(1, 5), a(1, 6), a(1, 7), a(1, 8), a(1, 9), a(1, 10)},
                {mx_i[1], a(2, 1), a(2, 2), a(2, 3), a(2, 4), a(2, 5), a(2, 6), a(2, 7), a(2, 8), a(2, 9), a(2, 10)},
                {mx_i[2], a(3, 1), a(3, 2), a(3, 3), a(3, 4), a(3, 5), a(3, 6), a(3, 7), a(3, 8), a(3, 9), a(3, 10)},
                {mx_i[3], a(4, 1), a(4, 2), a(4, 3), a(4, 4), a(4, 5), a(4, 6), a(4, 7), a(4, 8), a(4, 9), a(4, 10)},
                {mx_i[4], a(5, 1), a(5, 2), a(5, 3), a(5, 4), a(5, 5), a(5, 6), a(5, 7), a(5, 8), a(5, 9), a(5, 10)},
                {mx_i[5], a(6, 1), a(6, 2), a(6, 3), a(6, 4), a(6, 5), a(6, 6), a(6, 7), a(6, 8), a(6, 9), a(6, 10)},
                {mx_i[6], a(7, 1), a(7, 2), a(7, 3), a(7, 4), a(7, 5), a(7, 6), a(7, 7), a(7, 8), a(7, 9), a(7, 10)},
                {mx_i[7], a(8, 1), a(8, 2), a(8, 3), a(8, 4), a(8, 5), a(8, 6), a(8, 7), a(8, 8), a(8, 9), a(8, 10)},
                {mx_i[8], a(9, 1), a(9, 2), a(9, 3), a(9, 4), a(9, 5), a(9, 6), a(9, 7), a(9, 8), a(9, 9), a(9, 10)},
                {mx_i[9], a(10, 1), a(10, 2), a(10, 3), a(10, 4), a(10, 5), a(10, 6), a(10, 7), a(10, 8), a(10, 9), a(10, 10)
                }};

        double[] known = new double[]{my, findKnown(1, averageY), findKnown(2, averageY), findKnown(3, averageY), findKnown(4,  averageY), findKnown(5,  averageY), findKnown(6,  averageY), findKnown(7,  averageY),
                findKnown(8, averageY), findKnown(9, averageY), findKnown(10, averageY)};

        double[] b = solveMatrix.lsolve(unknown, known);
        List<Double> allB = new ArrayList<Double>();

        for(int i = 0; i < b.length; i++) {
            allB.add(b[i]);
//            System.out.println("b"+i+" = "+b[i]);
        }

        double[] allY = new double[N];
        for(int i = 0; i < allY.length; i++) {
            if(coef == 0) {
                allY[i] = allB.get(0) + allB.get(1) * matrixX[i][0] + allB.get(2) * matrixX[i][1] + +allB.get(3) * matrixX[i][2] + allB.get(4) * matrixX[i][3] +
                        +allB.get(5) * matrixX[i][4] + allB.get(6) * matrixX[i][5];
            }
            if(coef == 1) {
                allY[i] = allB.get(0) + allB.get(1) * matrixX[i][0] + allB.get(2) * matrixX[i][1] + +allB.get(3) * matrixX[i][2] + allB.get(4) * matrixX[i][3] +
                        +allB.get(5) * matrixX[i][4] + allB.get(6) * matrixX[i][5] + allB.get(7) * matrixX[i][6];
            }
            if(coef == 2) {
                allY[i] = allB.get(0) + allB.get(1) * matrixX[i][0] + allB.get(2) * matrixX[i][1] + +allB.get(3) * matrixX[i][2] + allB.get(4) * matrixX[i][3] +
                        +allB.get(5) * matrixX[i][4] + allB.get(6) * matrixX[i][5] + allB.get(7) * matrixX[i][6] + allB.get(8) * matrixX[i][7] + allB.get(9) * matrixX[i][8] + allB.get(10) * matrixX[i][9];
            }
        }

        //Перевірка
        if(bool){
            System.out.println("Перевірка:");
            for(int i = 0; i < averageY.length; i++){
                System.out.println("* y"+i+" = "+allY[i]+" ≈ "+averageY[i]);
            }
        }

        boolean homogeneity = false;
        while (!homogeneity){

            for(int i=0; i < N;i++){
                dispersionY.add(0.0);
            }

            for(int i=0; i < N; i++){
                double dispersionI = 0;
                for(int j=0; j < m; j++) {
                    dispersionI += Math.pow(matrixY[i][j] - averageY[i], 2);
                }
                dispersionY.add(dispersionI/(m-1));
            }
            int f1 = m -1;
            int f2 = N;
            f3 = f1 * f2;
            double q = 1 -p;
            double Gp = Collections.max(dispersionY) / sumList(dispersionY) ;
            double Gt = data.getTableCohren(f1, f2);

            if(Gt > Gp || m>=25){
                if(bool){
                    System.out.println("\nЗа критерієм Кохрена:\n* Дисперсія однорідна при рівні значимості : "+q);
                }
                homogeneity = true;
            }
            else {
                if(bool){
                    System.out.println("\nЗа критерієм Кохрена:\n* Дисперсія не однорідна при рівні значимості : "+q);
                }
                m+=1;
            }
            if(m==25){
                System.exit(0);
            }
        }

        dispersionB2 = sumList(dispersionY) / (N * N * m);
        List<Double> studList = studentTest(allB, 0);


        stY = new double[N];
        for(int i = 0; i < stY.length; i++) {
            if(coef==0) {
                stY[i] = studList.get(0) + studList.get(1) * matrixX[i][0] + studList.get(2) * matrixX[i][1] + studList.get(3) * matrixX[i][2] + studList.get(4) * matrixX[i][3] +
                        +studList.get(5) * matrixX[i][4] + studList.get(6) * matrixX[i][5];
            }
            if(coef==1) {
                stY[i] = studList.get(0) + studList.get(1) * matrixX[i][0] + studList.get(2) * matrixX[i][1] + studList.get(3) * matrixX[i][2] + studList.get(4) * matrixX[i][3] +
                        +studList.get(5) * matrixX[i][4] + studList.get(6) * matrixX[i][5] + studList.get(7) * matrixX[i][6];
            }
            if(coef==2) {
                stY[i] = studList.get(0) + studList.get(1) * matrixX[i][0] + studList.get(2) * matrixX[i][1] + studList.get(3) * matrixX[i][2] + studList.get(4) * matrixX[i][3] +
                        +studList.get(5) * matrixX[i][4] + studList.get(6) * matrixX[i][5] + studList.get(7) * matrixX[i][6] + studList.get(8) * matrixX[i][7] + studList.get(9) * matrixX[i][8] + studList.get(10) * matrixX[i][9];
            }
        }

        //Перевірка
        if(bool){
            System.out.println("\nЗа критерієм Стьюдента");
            for(int i = 0; i < averageY.length; i++){
                System.out.println("y"+(i)+" = "+stY[i]+" ≈ "+averageY[i]);
            }

            System.out.println("\nЗа критерієм Фішера");
        }
        int n0 = 0;
        for(int i=0; i < studList.size(); i++){
            if(studList.get(i) == 0.){
                n0++;
            }
        }
        d = 11 - n0;

        if(fisherTest()){
            if(bool){
                System.out.println("* Рівняння регресії адекватне стосовно оригіналу.");
            }
            printEquation(b);
            return true;
        } else {
            System.out.println("* Рівняння регресії неадекватне стосовно оригіналу.\nРівняння не буде...");
            return false;
        }

    }

    private static void printEquation(double[] b){
        if(coef == 1){
            System.out.printf("\nРівняння регресії з ефектом взаємодії:\ny = %+f%+f*X1%+f*X2%+f*X3%+f*X1*X2%+f*X1*X3" +
                    "%+f*X2*X3%+f*X1*X2*X3\n", b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]);
            coef = 2;
            nEf++;
        } else if(coef == 2){
            System.out.printf("\nРівняння регресії з квадратичними членами:\ny = %+f%+f*X1%+f*X2%+f*X3%+f*X1*X2%+f*X1*X3" +
                    "%+f*X2*X3%+f*X1*X2*X3%+f*X1*X1%+f*X2*X2%+f*X3*X3\n", b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],b[8],b[9],b[10] );
            coef = 0;
            nQt++;
        } else {
            System.out.printf("\nРівняння регресії лінійна форма:\ny = %+f%+f*X1%+f*X2%+f*X3\n", b[0], b[1], b[2], b[3], b[4]);
            coef = 1;
            nL++;
        }
    }

    public static double sum(double...values) {
        double result = 0;
        for (double value:values)
            result += value;
        return result;
    }

    public static double sumList(List<Double> list) {
        double sum = 0;
        for (double i: list) {
            sum += i;
        }
        return sum;
    }

    //зоряеі точки
    private static double[] countX(double l0,double l1,double l2){
        double x_1 = l0*deltaX1 + x01;
        double x_2 = l1*deltaX2 + x02;
        double x_3 = l2*deltaX3 + x03;

        return new double[]{x_1,x_2,x_3};
    }

    private static double[][] generateMatrix(int m, int n){
        double[][] matrix = new double[m][n];

        for(int i=0; i < m; i++){
            for(int j=0; j < n; j++){
                matrix[i][j] = minY + (maxY - minY)*r.nextDouble();
            }
        }
        return  matrix;
    }
    //k=0 -пошуксереднього по стовбцях к=1 - по рядках
    private static double[] getAverage(double[][] list, int k){
        double [] result;

        if(k == 0 && list.length != 0){
            result = new double[list[0].length];

            for (int i=0; i < list.length; i++){
                for (int j=0; j < list[i].length; j++){
                    result[j] += list[i][j];
                }
            }
            for(int i =0; i < result.length; i++){
                result[i] = result[i]/list.length;
            }

        }
        else {
            result = new double[list.length];

            double sumRow;
            for(int i=0; i < list.length; i++){
                sumRow = 0;
                for(int j=0; j < list[i].length; j++){
                    sumRow += list[i][j];
                }
                result[i] = sumRow/list[i].length;
            }
        }
        return result;
    };

    private static double a(int f, int s){
        double needA = 0;

        for(int i =0; i < N; i++){
            needA += matrixX[i][f-1]*matrixX[i][s-1]/N;
        }
        return needA;
    };

    private static double findKnown(int n, double[] average){
        double needA =0;

        for(int i =0; i < N; i++){
            needA += average[i]*matrixX[i][n-1]/N;
        }
        return needA;
    };

    private static void printM(double[][] list){
        for(int i = 0; i < list.length; i++){
            System.out.print("{");
            for (int j =0; j < list[i].length; j++){
                System.out.print(list[i][j]+ " ,");
            }
            System.out.print("}\n");
        }
    }
    //
    private static List<Double> studentTest(List<Double> allB, int nx){
        nx = nx == 0 ? 10 : nx;
        double dispersionB = Math.sqrt(dispersionB2);

        for(int c=0; c<nx; c++){
            double tPractice =0;
            double tTheoretical = data.getTableStudent(f3, p);

            for(int r=0; r < N; r++){
                if(c ==  0){
                    tPractice += averageY[r]/N;
                }
                else{
                    tPractice += averageY[r]*normMatrix[r][c-1];
                }
            }
            if(Math.abs(tPractice/dispersionB) < tTheoretical){
                allB.set(c, 0.);
            }
        }
        return allB;
    }

    private static boolean fisherTest(){
        double dispersionAb = 0;
        int f4 = N - d;

        for(int i =0; i < averageY.length; i++){
            dispersionAb += (m*(averageY[i] - stY[i])) / (N-d);
        }
        double practiceF = Math.abs(dispersionAb / dispersionB2);
        double theoreticalF = data.getFisherValue(f3,f4,q);
        try {
            theoreticalF  = Float.parseFloat(System.console().readLine());
        }
        catch (Exception ex){
        }
        if(bool){
            if(practiceF < theoreticalF){
                System.out.println("* F_пр < F_тр");
            } else {
                System.out.println("* F_пр > F_тр");
            }
        }
        return practiceF < theoreticalF;
    }

}

class GaussianElimination {
    private final double EPSILON = 1e-10;

    // Gaussian elimination with partial pivoting
    public double[] lsolve(double[][] A, double[] b) {
        int n = b.length;

        for (int p = 0; p < n; p++) {

            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < n; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }
            double[] temp = A[p]; A[p] = A[max]; A[max] = temp;
            double   t    = b[p]; b[p] = b[max]; b[max] = t;

            // singular or nearly singular
            if (Math.abs(A[p][p]) <= EPSILON) {
                throw new ArithmeticException("Matrix is singular or nearly singular");
            }

            // pivot within A and b
            for (int i = p + 1; i < n; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < n; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }

        // back substitution
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }
}

class Data {

    public double[][] cohrenCriterium = {
            {9985, 9750, 9392, 9057, 8772, 8534},
            {9669, 8709, 7977, 7457, 7071, 6771},
            {9065, 7679, 6841, 6282, 5892, 5598},
            {8412, 6838, 5981, 5440, 5063, 4783},
            {7808, 6161, 5321, 4803, 4447, 4184},
            {7271, 5612, 4800, 4307, 3974,3726},
            {6798, 5157, 4377, 3910, 3595, 3362},
            {6385,4775, 4027, 3584, 3286, 3067},
            {6020, 4450, 3733, 3311, 3029, 2823}
    };

    // 0.999 0.9 0.99 0.90
    public double[][] studentTable = {
            {6.3137515148,	12.7062047364,	63.6567411629,	636.619249432},
            {2.91998558036,	4.30265272991,	9.92484320092,	31.599054577},
            {2.3533634348,	3.18244630528,	5.84090929976,	12.9239786366},
            {2.13184678134,	2.7764451052,	4.60409487142,	8.61030158138},
            {2.01504837267,	2.57058183661,	4.03214298356,	6.86882663987},
            {1.94318028039,	2.44691184879,	3.70742802132,	5.95881617993},
            {1.89457860506,	2.36462425101,	3.49948329735,	5.40788252098},
            {1.85954803752,	2.30600413503,	3.35538733133,	5.04130543339},
            {1.83311293265,	2.26215716274,	3.24983554402,	4.78091258593},
            {1.81246112281,	2.22813885196,	3.16927266718,	4.5868938587},
            {1.7958848187,	2.20098516008,	3.10580651322,	4.43697933823},
            {1.78228755565,	2.17881282966,	3.05453958834,	4.31779128361},
            {1.77093339599,	2.16036865646,	3.01227583821,	4.22083172771},
            {1.76131013577,	2.14478668792,	2.97684273411,	4.14045411274},
            {1.75305035569,	2.13144954556,	2.94671288334,	4.0727651959},
            {1.74588367628,	2.11990529922,	2.92078162235,	4.0149963326},
            {1.73960672608,	2.10981557783,	2.89823051963,	3.96512626361},
            {1.73406360662,	2.10092204024,	2.87844047271,	3.92164582001},
            {1.72913281152,	2.09302405441,	2.86093460645,	3.88340584948},
            {1.72471824292,	2.08596344727,	2.84533970978,	3.84951627298},
            {1.72074290281,	2.07961384473,	2.83135955802,	3.81927716303},
            {1.71714437438,	2.0738730679,	2.8187560606,	3.79213067089},
            {1.71387152775,	2.06865761042,	2.80733568377,	3.76762680377},
            {1.71088207991,	2.06389856163,	2.79693950477,	3.74539861893},
            {1.70814076125,	2.05953855275,	2.78743581368,	3.72514394948},
            {1.70561791976,	2.05552943864,	2.77871453333,	3.70661174331},
            {1.70328844572,	2.05183051648,	2.77068295712,	3.68959171334},
            {1.70113093427,	2.0484071418,	2.76326245546,	3.67390640062},
            {1.69912702653,	2.04522964213,	2.75638590367,	3.6594050194},
            {1.69726089436,	2.0422724563,	2.74999565357,	3.645958635},
            {1.68385101139,	2.021075383,	2.70445926743,	3.55096576086},
            {1.67064886465,  2.00029782106,	2.66028303115,	3.4602004692},
            {1.65765089935,	1.97993040505,	2.61742114477,	3.37345376507},
            {1000000.0,	1.64485515072,	1.95996635682,	2.57583422011,	3.29053646126}
    };

    //30 45 60
    public double[] fisherTable =
            {2.16 , 2.06, 2.036};

    public double getTableCohren(int f1, int f2){
        double[] n = new double[]{1,2,3,4,5,6,7,8,9,10,16,36,144};

        for(int i=0; i < n.length; i++){
            if(f2 == n[i]){
                return cohrenCriterium[f1][i];
            }
            if(f2 < n[i] || i>1){
                double c1 = cohrenCriterium[f1][i-1];
                double c2 = cohrenCriterium[f1][i];
                return c1+(c2-c1)/2;
            }
        }
        return 0.;
    }

    public double getTableStudent(int f3, double q){
        double[] allQ = new double[]{0.90, 0.95, 0.99, 0.999};

        for(int i=0; i < allQ.length; i++) {
            if (q == allQ[i]) {
                if( f3 > 30 && f3 <= 40){
                    return studentTable[30][i];
                }
                if(f3 > 40 && f3 <= 60){
                    return studentTable[31][i];
                }
                if(f3 > 60 && f3 <= 120){
                    return studentTable[32][i];
                }
                if(f3 > 120){
                    return studentTable[33][i];
                }
                return studentTable[f3][i];
            }
        }
        return 0.;
    }
    public double getFisherValue(int f3, int f4, double q){
        if(f3==30)
            return fisherTable[0];
        if(f3 == 45)
            return fisherTable[1];
        if(f3 == 60)
            return fisherTable[2];
        return 164.;
    }
}

