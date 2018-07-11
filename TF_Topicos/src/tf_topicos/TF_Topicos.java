package tf_topicos;

import java.util.*;
import java.io.*;
import weka.core.*;
import weka.core.converters.ConverterUtils;

public class TF_Topicos {

    private boolean[] marked;
    private int[] edgeTo;
    private boolean[] onStack;
    private Stack<Integer> Cycle;

    public TF_Topicos(int[][] G) {
        marked = new boolean[G.length];
        onStack = new boolean[G.length];
        edgeTo = new int[G.length];
        for (int v = 0; v < G.length; v++) {
            if (!marked[v] && Cycle == null) {
                dfs(G, v);
            }
        }
    }

    public int[] obtenerVecinos(int[][] G, int v) {
        int i;
        List<Integer> hijos = new ArrayList<Integer>();
        for (i = 0; i < G.length; i++) {
            if (G[i][v] != 0) {
                hijos.add(i);
            }
        }
        int[] nums = new int[hijos.size()];
        for (i = 0; i < hijos.size(); i++) {
            nums[i] = hijos.get(i);
        }
        return nums;
    }

    public void dfs(int[][] G, int v) {
        onStack[v] = true;
        marked[v] = true;
        int[] vecino = obtenerVecinos(G, v);
        for (Integer w : vecino) {
            if (Cycle != null) {
                return;
            } else if (!marked[w]) {
                edgeTo[w] = v;
                dfs(G, w);
            } else if (onStack[w]) {
                Cycle = new Stack<Integer>();
                for (int x = v; x != w; x = edgeTo[x]) {
                    Cycle.push(x);
                }
                Cycle.push(w);
                Cycle.push(v);
            }
        }
        onStack[v] = false;
    }

    public boolean hasCycle() {
        if (Cycle == null) {
            return false;
        }
        return true;
    }

    public static int[] HipotesisMatrix(int n, int numOfBits) {
        int[] binary = new int[numOfBits];
        for (int i = 0; i < numOfBits; ++i, n /= 2) {
            switch (n % 2) {
                case 0:
                    binary[i] = 0;
                    break;
                case 1:
                    binary[i] = 1;
                    break;
            }
        }
        return binary;
    }

    public static List<Integer> getPadres(int i, int[][] grafo) {
        List<Integer> padres = new ArrayList<>();
        for (int k = 0; k < grafo.length; k++) {
            if (grafo[i][k] == 1) {
                padres.add(k);
            }
        }
        return padres;
    }

    public static List<List> getInferencia(int[][] matrix) {
        //System.out.print("Creando inferencia");
        List<List> myList = new ArrayList<List>();
        for (int i = 0; i < matrix.length; i++) {
            List<Integer> temp = new ArrayList<>();
            temp.add(i);
            for (int j = 0; j < matrix.length; j++) {
                if (matrix[j][i] == 1) {
                    temp.add(j);
                }

            }
            myList.add(temp);
        }
        return myList;
    }

    public static double[] maximo(double[] inf) {
        double max = 0.0;
        int maxid = 0;
        double[] prob = new double[2];
        for (int i = 0; i < inf.length; i++) {
            if (inf[i] > max) {
                maxid = i;
                max = inf[i];
                prob[0] = max;
                prob[1] = maxid;
            }
        }

        return prob;
    }

    public static int CalcularIndice(int[] card, int[] vars, int[] valinf, int l) {

        int indice = 0;
        for (int i = 0; i < vars.length; i++) {
            int cardAnt = 1;

            for (int j = 0; j < i; j++) {
                cardAnt *= card[vars[j]];

            }
            if (valinf[vars[i]] != -1) {
                cardAnt *= valinf[vars[i]];
            } else {
                cardAnt *= l;
            }

            indice += cardAnt;
        }

        return indice;

    }

    public static int[] union(int i, List<Integer> padres) {
        int N = padres.size() + 1;
        int[] arr = new int[N];
        int k = 0;
        arr[k++] = i;
        for (Integer x : padres) {
            arr[k++] = x;
        }
        return arr;
    }

    public static int TablaConteoMarginal(int var, int[][] data, int val) {
        int prob = 0;
        if (data.length > 0) {
            for (int i = 0; i < data.length; i++) {
                if (data[i][var] == val) {
                    prob++;
                }
            }
        }
        return prob;
    }

    public static int TablaConteoConjunta(int[] vars, int[][] data, int[] vals) {
        int contador = 0;
        if (data.length > 0) {
            boolean boleano = false;
            double cardmul = 1.0;
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < vars.length; j++) {
                    if (data[i][vars[j]] != vals[j]) {
                        boleano = false;
                        break;
                    } else {
                        boleano = true;
                    }
                }
                if (boleano) {
                    contador++;
                }
            }
        }
        return contador;
    }

    public static List<int[]> FactorConjunta(int N, int i) {
        List<int[]> combs = new ArrayList<int[]>();
        int r = 0;
        int[] comb = new int[i];
        int index = 0;
        while (r >= 0) {
            if (index <= (N + (r - i))) {
                comb[r] = index;
                if (r == i - 1) {
                    int[] aux = comb.clone();
                    combs.add(aux);
                    index++;
                } else {
                    index = comb[r] + 1;
                    r++;
                }
            } else {
                r--;
                if (r > 0) {
                    index = comb[r] + 1;
                } else {
                    index = comb[0] + 1;
                }
            }
        }
        return combs;
    }

    public static List<int[]> FactorCondicional(int N, int i) {
        List<int[]> combs = new ArrayList<int[]>();
        int r = 0;
        int[] comb = new int[i];
        int index = 0;
        while (r >= 0) {
            if (index <= N - 1) {
                comb[r] = index;
                if (r == i - 1) {
                    boolean t = true;
                    HashSet<Integer> hs = new HashSet<>();
                    for (int num : comb) {
                        if (hs.contains(num)) {
                            t = false;
                            break;
                        }
                        hs.add(num);
                    }
                    if (t) {
                        for (int j = 0; j < combs.size(); j++) {
                            int[] aux2 = combs.get(j);
                            if (aux2.length != comb.length || aux2[0] != comb[0]) {
                                continue;
                            }
                            HashSet<Integer> hs2 = new HashSet<>();
                            for (int num : aux2) {
                                hs2.add(num);
                            }
                            if (hs.equals(hs2)) {
                                t = false;
                                break;
                            }
                        }

                    }
                    if (t) {
                        int[] aux = comb.clone();
                        combs.add(aux);
                        index++;
                    } else {
                        index++;
                    }
                } else {
                    index = 0;
                    r++;
                }
            } else {
                r--;
                if (r > 0) {
                    index = comb[r] + 1;
                } else {
                    index = comb[0] + 1;
                }
            }
        }
        return combs;
    }

    public static double CalcularProbMarginalDirichlet(int var, int[][] data, int val, int[] card, double alpha) {
        double prob = 0.0;
        if (data.length > 0) {
            for (int i = 0; i < data.length; i++) {
                if (data[i][var] == val) {
                    prob++;
                }
            }
            prob = (prob + alpha) / (data.length + card[var] * alpha);
        }
        return prob;
    }

    public static double CalcularProbConjuntaDirichlet(int[] vars, int[][] data, int[] vals, int[] card, double alpha, boolean condicional) {
        double proba = 0.0;
        if (data.length > 0) {
            double contador = 0;
            boolean boleano = false;
            double cardmul = 1.0;
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < vars.length; j++) {
                    if (data[i][vars[j]] != vals[j]) {
                        boleano = false;
                        break;
                    } else {
                        boleano = true;
                    }
                }
                if (boleano) {
                    contador++;
                }
            }

            for (int i = 0; i < vals.length; i++) {
                cardmul = cardmul * card[vars[i]];
            }

            if (condicional && contador == 0) {
                contador = (contador + alpha) / ((data.length) + (cardmul * alpha));
                return contador;
            }

            proba = (contador + alpha) / ((data.length) + (cardmul * alpha));
        }
        return proba;
    }

    public static double[] CalcularProbCondicionalDirichlet(int[] vars, int[][] data, int[][] vals, int[] card) {
        double proba = 0.0;
        double probb = 0.0;
        double[] probs = new double[card[vars[0]]];
        int[] aux = Arrays.copyOfRange(vars, 1, vars.length);
        int[] aux2;
        for (int i = 0; i < card[vars[0]]; i++) {
            aux2 = Arrays.copyOfRange(vals[i], 1, vals[i].length);
            proba = CalcularProbConjuntaDirichlet(vars, data, vals[i], card, 1.0, true);
            if (proba == 0) {
                probs[i] = proba;
                break;
            }
            if (aux.length > 1) {
                probb = CalcularProbConjuntaDirichlet(aux, data, aux2, card, 1.0, true);
            } else {
                probb = CalcularProbMarginalDirichlet(aux[0], data, aux2[0], card, 1.0);
            }
            probs[i] = proba / probb;
        }
        double probsum = 0;
        double total = 1.0;
        int auxcont = 0;
        for (int i = 0; i < probs.length; i++) {
            if (probs[i] == 0) {
                auxcont++;
            }
            probsum += probs[i];
        }
        if (probsum != 1.0) {
            total -= probsum;
            total /= (auxcont);
            for (int i = 0; i < probs.length; i++) {
                if (probs[i] == 0) {
                    probs[i] = total;
                }
            }
        }
        probsum = 0;
        for (int i = 0; i < probs.length; i++) {
            probsum += probs[i];
        }
        //System.out.println(probsum);

        return probs;
    }

    public static double CalcularProbMarginalEntropia(int var, int[][] data, int val, int[] card, double alpha) {
        double prob = 0.0;
        double cont = 0.0;
        if (data.length > 0) {
            for (int i = 0; i < data.length; i++) {
                if (data[i][var] == val) {
                    prob++;
                }
            }
            cont = TablaConteoMarginal(var, data, val);
            prob = (prob + alpha) / (data.length + card[var] * alpha);
            prob = cont * Math.log10(prob);
        }
        return prob;
    }

    public static double CalcularMarginalEntropia(int[] var, int[][] data, int[] val, int[] card, double alpha) {

        double prob = 0.0;
        if (data.length > 0) {
            double contador = 0;
            boolean boleano = false;
            double cardmul = 1.0;
            for (int i = 0; i < data.length; i++) {
                for (int j = 1; j < var.length; j++) {
                    if (data[i][var[j]] != val[j]) {
                        boleano = false;
                        break;
                    } else {
                        boleano = true;
                    }
                }

                if (boleano) {
                    contador++;
                }
            }

            for (int i = 1; i < val.length; i++) {
                cardmul = cardmul * card[var[i]];
            }

            if (contador == 0) {
                contador = (contador + alpha) / ((data.length) + (cardmul * alpha));
                return contador;
            }

            prob = (contador + alpha) / ((data.length) + (cardmul * alpha));;
        }
        return prob;

    }

    public static int getIndex(int[] data, int[] card) {
        int size = card.length;
        int sizeL = 1;
        for (int i = 0; i < card.length; i++) {
            sizeL *= card[i];
        }
        int position = 0;
        int half = sizeL / card[size - 1];
        for (int i = size - 1; i > 0; i--) {
            position += half * data[i];

            half = half / card[i - 1];

        }
        position += data[0];
        return position;
    }

    public static double[] CalcularProbCondicionalEntropia(int[] vars, int[][] data, int[][] vals, int[] card) {
        double proba = 0.0;
        double probb = 0.0;
        int conteo = 0;
        double[] probs = new double[card[vars[0]]];
        int[] aux = Arrays.copyOfRange(vars, 0, vars.length);
        int[] aux2;
        for (int i = 0; i < card[vars[0]]; i++) {
            aux2 = Arrays.copyOfRange(vals[i], 0, vals[i].length);

            //proba = CalcularProbConjuntaDirichlet(vars, data, vals[i], card, 1.0, true);
            //probb = CalcularProbMarginalDirichlet(aux[0], data, aux2[0], card, 1.0);
            if (aux.length > 1) {

                proba = CalcularProbConjuntaDirichlet(vars, data, vals[i], card, 1.0, true);
                probb = CalcularMarginalEntropia(aux, data, aux2, card, 1.0);
                //probb = CalcularProbMarginalDirichlet(aux[0], data, aux2[0], card, 1.0);
                probs[i] = proba / probb;
            } else {

                probb = CalcularProbMarginalDirichlet(aux[0], data, aux2[0], card, 1.0);
                probs[i] = probb;
            }

            /*if (proba == 0) {
                probs[i] = proba;
                break;
            }
            if (aux.length > 1) {
                probb = CalcularProbConjuntaDirichlet(aux, data, aux2, card, 1.0, true);
                probs[i] = proba / probb;
            } else {
                probb = CalcularProbMarginalDirichlet(aux[0], data, aux2[0], card, 1.0);
                probs[i] = probb;
            }*/
            //conteo = TablaConteoConjunta(vars, data, vals[i]);
            //probs[i] = conteo * Math.log10(probs[i]);
        }

        double total3 = 0.0;
        double suma = 0.0;
        for (int h = 0; h < card[vars[0]]; h++) {

            total3 = probs[h];
            suma = suma + total3;

        }
        for (int i = 0; i < card[vars[0]]; i++) {
            probs[i] = probs[i] / suma;
        }

        for (int i = 0; i < probs.length; i++) {
            conteo = TablaConteoConjunta(vars, data, vals[i]);
            probs[i] = conteo * Math.log10(probs[i]);
        }

        return probs;
    }

    public static double[] CalcularFactor(int iFactor, double[] probs, int[] vars, int[] card, int[][] data, String[] varnames, String[][] varvals, boolean condicional) throws FileNotFoundException {
        int N = probs.length;
        int m = vars.length;
        int[] vals = new int[m];
        int cardAnt = 1;
        int iVar;

        if (!condicional) {
            for (int i = 0; i < N; i++) {
                cardAnt = 1;
                for (int j = 0; j < m; j++) {
                    iVar = vars[j];
                    vals[j] = (int) (Math.floor(i / cardAnt) % card[iVar]);
                    cardAnt *= card[iVar];
                }
                double prob = 0.0;
                if (m > 1) {
                    prob = CalcularProbConjuntaDirichlet(vars, data, vals, card, 1.0, false);
                } else {
                    prob = CalcularProbMarginalDirichlet(vars[0], data, vals[0], card, 1.0);
                }
                probs[i] = prob;
            }
        } else {
            int[][] auxvals = new int[card[vars[0]]][m];
            for (int i = 0; i < N; i += card[vars[0]]) {
                for (int k = 0; k < card[vars[0]]; k++) {
                    cardAnt = 1;
                    for (int j = 0; j < m; j++) {
                        iVar = vars[j];
                        vals[j] = (int) (Math.floor((i + k) / cardAnt) % card[iVar]);
                        cardAnt *= card[iVar];
                    }
                    auxvals[k] = vals.clone();
                }
                double[] prob;
                prob = CalcularProbCondicionalDirichlet(vars, data, auxvals, card);

            }
        }
        return probs;
    }

    public static double[] ShowCalcularFactor(int iFactor, double[] probs, int[] vars, int[] card, int[][] data, String[] varnames, String[][] varvals, boolean condicional) throws FileNotFoundException {
        int N = probs.length;
        int m = vars.length;
        int[] vals = new int[m];
        int cardAnt = 1;
        int iVar;

        System.out.println("Factor " + iFactor);
        int w;
        System.out.print("i\t");
        for (w = 0; w < vars.length - 1; w++) {
            System.out.print(varnames[vars[w]] + "\t");
        }
        System.out.println(varnames[vars[w]]);
        if (!condicional) {
            for (int i = 0; i < N; i++) {
                cardAnt = 1;
                for (int j = 0; j < m; j++) {
                    iVar = vars[j];
                    vals[j] = (int) (Math.floor(i / cardAnt) % card[iVar]);
                    cardAnt *= card[iVar];
                }
                double prob = 0.0;
                if (m > 1) {
                    prob = CalcularProbConjuntaDirichlet(vars, data, vals, card, 1.0, false);
                } else {
                    prob = CalcularProbMarginalDirichlet(vars[0], data, vals[0], card, 1.0);
                }
                probs[i] = prob;
                System.out.print(i + "\t");
                for (w = 0; w < vals.length - 1; w++) {
                    System.out.print(varvals[vars[w]][vals[w]] + "\t");
                }
                System.out.println(varvals[vars[w]][vals[w]] + "\t" + prob);
            }
        } else {
            int[][] auxvals = new int[card[vars[0]]][m];
            for (int i = 0; i < N; i += card[vars[0]]) {
                for (int k = 0; k < card[vars[0]]; k++) {
                    cardAnt = 1;
                    for (int j = 0; j < m; j++) {
                        iVar = vars[j];
                        vals[j] = (int) (Math.floor((i + k) / cardAnt) % card[iVar]);
                        cardAnt *= card[iVar];
                    }
                    auxvals[k] = vals.clone();
                }
                double[] prob;
                prob = CalcularProbCondicionalDirichlet(vars, data, auxvals, card);
                for (int l = 0; l < card[vars[0]]; l++) {

                    double total = 0.0;
                    double suma = 0.0;

                    for (int h = 0; h < card[vars[0]]; h++) {

                        total = prob[h];
                        suma = suma + total;

                    }

                    probs[i + l] = prob[l] / suma;

                    System.out.print((i + l) + "\t");

                    for (w = 0; w < vals.length - 1; w++) {
                        System.out.print(varvals[vars[w]][auxvals[l][w]] + "\t");
                    }
                    System.out.println(varvals[vars[w]][auxvals[l][w]] + "\t" + prob[l] / suma);
                }
            }
        }
        return probs;
    }

    public static double[] CalcularProbCondicionalResul(int[] vars, int[][] data, int[][] vals, int[] card) {
        double proba = 0.0;
        double probb = 0.0;
        int conteo = 0;
        double[] probs = new double[card[vars[0]]];
        int[] aux = Arrays.copyOfRange(vars, 0, vars.length);
        int[] aux2;
        for (int i = 0; i < card[vars[0]]; i++) {
            aux2 = Arrays.copyOfRange(vals[i], 0, vals[i].length);

            //proba = CalcularProbConjuntaDirichlet(vars, data, vals[i], card, 1.0, true);
            //probb = CalcularProbMarginalDirichlet(aux[0], data, aux2[0], card, 1.0);
            if (aux.length > 1) {

                proba = CalcularProbConjuntaDirichlet(vars, data, vals[i], card, 1.0, true);
                probb = CalcularMarginalEntropia(aux, data, aux2, card, 1.0);
                //probb = CalcularProbMarginalDirichlet(aux[0], data, aux2[0], card, 1.0);
                probs[i] = proba / probb;
            } else {

                probb = CalcularProbMarginalDirichlet(aux[0], data, aux2[0], card, 1.0);
                probs[i] = probb;
            }

            /*if (proba == 0) {
                probs[i] = proba;
                break;
            }
            if (aux.length > 1) {
                probb = CalcularProbConjuntaDirichlet(aux, data, aux2, card, 1.0, true);
                probs[i] = proba / probb;
            } else {
                probb = CalcularProbMarginalDirichlet(aux[0], data, aux2[0], card, 1.0);
                probs[i] = probb;
            }*/
            //conteo = TablaConteoConjunta(vars, data, vals[i]);
            //probs[i] = conteo * Math.log10(probs[i]);
        }

        double total3 = 0.0;
        double suma = 0.0;
        for (int h = 0; h < card[vars[0]]; h++) {

            total3 = probs[h];
            suma = suma + total3;

        }
        for (int i = 0; i < card[vars[0]]; i++) {
            probs[i] = probs[i] / suma;
        }

        /*for (int i = 0; i < probs.length; i++) {
            conteo = TablaConteoConjunta(vars, data, vals[i]);
            probs[i] = conteo * Math.log10(probs[i]);
        }*/
        return probs;
    }

    public static double[] CalcularResul(int iFactor, double[] probs, int[] vars, int[] card, int[][] data, boolean condicional) {
        int N = probs.length;
        int m = vars.length;
        List<double[]> factores = new ArrayList<>();
        int[] vals = new int[m];
        int cardAnt = 1;
        int iVar;
        if (!condicional) {
            for (int i = 0; i < N; i++) {
                cardAnt = 1;
                for (int j = 0; j < m; j++) {
                    iVar = vars[j];
                    vals[j] = (int) (Math.floor(i / cardAnt) % card[iVar]);
                    cardAnt *= card[iVar];
                }
                double prob = 0.0;
                prob = CalcularProbMarginalEntropia(vars[0], data, vals[0], card, 1.0);
                probs[i] = prob;
            }
        } else {
            int[][] auxvals = new int[card[vars[0]]][m];
            int variable = 0;
            for (int i = 0; i < N; i += card[vars[0]]) {
                for (int k = 0; k < card[vars[0]]; k++) {
                    variable++;
                    cardAnt = 1;
                    for (int j = 0; j < m; j++) {
                        iVar = vars[j];
                        vals[j] = (int) (Math.floor((i + k) / cardAnt) % card[iVar]);
                        cardAnt *= card[iVar];
                    }
                    auxvals[k] = vals.clone();
                }
                double[] prob;
                prob = CalcularProbCondicionalResul(vars, data, auxvals, card);
                //factores.add(prob);*/

                for (int l = 0; l < card[vars[0]]; l++) {
                    probs[i + l] = prob[l];
                }

            }

        }

        return probs;
    }

    public static double[] CalcularEntropia(double[] probs, int[] vars, int[] card, int[][] data, boolean condicional) {
        int N = probs.length;
        int m = vars.length;
        int[] vals = new int[m];
        int cardAnt = 1;
        int iVar;
        if (!condicional) {
            for (int i = 0; i < N; i++) {
                cardAnt = 1;
                for (int j = 0; j < m; j++) {
                    iVar = vars[j];
                    vals[j] = (int) (Math.floor(i / cardAnt) % card[iVar]);
                    cardAnt *= card[iVar];
                }
                double prob = 0.0;
                prob = CalcularProbMarginalEntropia(vars[0], data, vals[0], card, 1.0);
                probs[i] = prob;
            }
        } else {
            int[][] auxvals = new int[card[vars[0]]][m];
            for (int i = 0; i < N; i += card[vars[0]]) {
                for (int k = 0; k < card[vars[0]]; k++) {
                    cardAnt = 1;
                    for (int j = 0; j < m; j++) {
                        iVar = vars[j];
                        vals[j] = (int) (Math.floor((i + k) / cardAnt) % card[iVar]);
                        cardAnt *= card[iVar];
                    }
                    auxvals[k] = vals.clone();
                }
                double[] prob;
                prob = CalcularProbCondicionalEntropia(vars, data, auxvals, card);

                for (int l = 0; l < card[vars[0]]; l++) {
                    probs[l + i] = prob[l];
                }

                /*for (int l = 0; l < card[vars[0]]; l++) {

                    double total = 0.0;
                    double suma = 0.0;

                    for (int h = 0; h < card[vars[0]]; h++) {

                        total = prob[h];
                        suma = suma + total;

                    }

                    probs[i + l] = prob[l] / suma;
                }*/
            }
        }
        return probs;
    }

    public static double[] CalcularConteo(int iFactor, double[] probs, int[] vars, int[] card, int[][] data, String[] varnames, String[][] varvals, boolean condicional) throws FileNotFoundException {
        int N = probs.length;
        int m = vars.length;
        int[] vals = new int[m];
        int cardAnt = 1;
        int iVar;
        System.out.println("Factor " + iFactor);
        int w;
        System.out.print("i\t");
        for (w = 0; w < vars.length - 1; w++) {
            System.out.print(varnames[vars[w]] + "\t");
        }
        System.out.println(varnames[vars[w]]);
        if (!condicional) {
            for (int i = 0; i < N; i++) {
                cardAnt = 1;
                for (int j = 0; j < m; j++) {
                    iVar = vars[j];
                    vals[j] = (int) (Math.floor(i / cardAnt) % card[iVar]);
                    cardAnt *= card[iVar];
                }
                int prob = 0;
                if (m > 1) {
                    prob = TablaConteoConjunta(vars, data, vals);
                } else {
                    prob = TablaConteoMarginal(vars[0], data, vals[0]);
                }
                probs[i] = prob;
                System.out.print(i + "\t");
                for (w = 0; w < vals.length - 1; w++) {
                    System.out.print(varvals[vars[w]][vals[w]] + "\t");
                }
                System.out.println(varvals[vars[w]][vals[w]] + "\t" + prob);
            }
        }
        return probs;
    }

    public static int GetAkaike(List<int[]> vars, int card[]) {
        int ka = 0;
        int cardp;
        int j = 0;
        for (int i = 0; i < vars.size(); i++) {
            cardp = 1;

            for (j = 1; j < vars.get(i).length; j++) {
                cardp *= card[vars.get(i)[j]];
            }
            ka += (card[i] - 1) * cardp;
        }

        return ka;
    }

    public static double GetEntropia(int[][] result, int[] card, int[][] datatrain, int N) {
        double e = 0.0;
        for (int i = 0; i < N; i++) {
            List<Integer> padres1 = getPadres(i, result);
            int[] vars1 = union(i, padres1);
            int tamFactor1 = 1;
            for (int l = 0; l < vars1.length; l++) {
                tamFactor1 *= card[vars1[l]];
            }
            double[] probs = new double[tamFactor1];
            probs = CalcularEntropia(probs, vars1, card, datatrain, true);
            if (vars1.length > 1) {
                for (int b = 0; b < probs.length; b++) {
                    e += probs[b];
                }
            } else {
                for (int b = 0; b < probs.length; b++) {
                    e += probs[b];
                }
            }
        }
        return e;
    }
    
    public static List<int[][]> GetHC(int[][] result, int[][] result2, int[][] result3, int[][] result4, TF_Topicos dc4, TF_Topicos dc2, TF_Topicos dc3) {
        List<int[][]> ArreglosHC = new ArrayList<int[][]>();

        do {

            for (int ip = 0; ip < result.length; ip++) {
                for (int ip2 = 0; ip2 < result.length; ip2++) {
                    result2[ip][ip2] = result[ip][ip2];
                }
            }

            int conta = 0;
            for (int i2 = 0; i2 < result.length; i2++) {
                for (int i3 = 0; i3 < result.length; i3++) {
                    if (result[i2][i3] == 0) {
                        conta++;
                    }
                }
            }

            if (conta <= 18) {
                int numero = (int) (Math.random() * 6);
                int numero2 = (int) (Math.random() * 6);
                while (result[numero][numero2] == 1 && numero == numero2) {
                    numero = (int) (Math.random() * 6);
                    numero2 = (int) (Math.random() * 6);
                }

                result2[numero][numero2] = 1;
                dc2 = new TF_Topicos(result2);
            } else {
                break;
            }

        } while (dc2.hasCycle());

        do {

            for (int ip = 0; ip < result.length; ip++) {
                for (int ip2 = 0; ip2 < result.length; ip2++) {
                    result3[ip][ip2] = result[ip][ip2];
                }
            }
            int conta = 0;
            for (int i2 = 0; i2 < result.length; i2++) {
                for (int i3 = 0; i3 < result.length; i3++) {
                    if (result[i2][i3] == 1) {
                        conta++;
                    }
                }
            }

            if (conta <= 18) {
                int numero = (int) (Math.random() * 6);
                int numero2 = (int) (Math.random() * 6);
                int cuenta = 0;
                while (result[numero][numero2] == 0 && numero == numero2) {

                    numero = (int) (Math.random() * 6);
                    numero2 = (int) (Math.random() * 6);
                }

                result3[numero][numero2] = 0;
                dc3 = new TF_Topicos(result3);
            } else {
                break;
            }

        } while (dc3.hasCycle());

        do {

            for (int ip = 0; ip < result.length; ip++) {
                for (int ip2 = 0; ip2 < result.length; ip2++) {
                    result4[ip][ip2] = result[ip][ip2];
                }
            }
            int conta = 0;
            for (int i2 = 0; i2 < result.length; i2++) {
                for (int i3 = 0; i3 < result.length; i3++) {
                    if (result[i2][i3] == 1) {
                        conta++;
                    }
                }
            }

            if (conta <= 18) {
                int numero = (int) (Math.random() * 6);
                int numero2 = (int) (Math.random() * 6);
                int cuenta = 0;
                while (result[numero][numero2] == 0 && numero == numero2) {

                    numero = (int) (Math.random() * 6);
                    numero2 = (int) (Math.random() * 6);
                }

                result4[numero2][numero] = 1;
                result4[numero][numero2] = 0;
                dc4 = new TF_Topicos(result4);
            } else {
                break;
            }

        } while (dc4.hasCycle());

        ArreglosHC.add(result);
        ArreglosHC.add(result2);
        ArreglosHC.add(result3);
        ArreglosHC.add(result4);

        return ArreglosHC;
    }

    public static List<String[]> Leer(File File) {
        List<String[]> lineas = new ArrayList<String[]>();
        try {
            FileReader fr = new FileReader(File);
            BufferedReader br = new BufferedReader(fr);
            String str;

            while ((str = br.readLine()) != null) {
                String[] lineItems = str.split(",");
                lineas.add(lineItems);

            }
            br.close();
        } catch (IOException e) {
            System.out.println("File not found");
        }
        return lineas;
    }

    public static List<int[]> GetVars(int N, int[][] G) {
        List<int[]> v = new ArrayList<>();
        int[] vars;
        for (int j = 0; j < N; j++) {
            List<Integer> padres = getPadres(j, G);
            vars = union(j, padres);
            v.add(vars);
        }
        return v;
    }

    public static List<double[]> GetProbabilidades(int[][] result, int[][] datatrain, int N, int[] card) {
        List<double[]> Probabilidades = new ArrayList<>();
        for (int q = 0; q < N; q++) {
            List<Integer> padres = getPadres(q, result);
            int[] vars = union(q, padres);
            int tamFactor = 1;
            for (int l = 0; l < vars.length; l++) {
                tamFactor *= card[vars[l]];
            }
            double[] probs = new double[tamFactor];
            double[] cont = new double[tamFactor];
            double[] probsNorm = new double[tamFactor];
            double res = 0.0;
            probs = CalcularResul(q, probs, vars, card, datatrain, true);
            if (vars.length > 1) {
                for (int b = 0; b < probs.length; b++) {
                    probsNorm[b] = probs[b];
                }
                Probabilidades.add(probsNorm);
            } else {
                Probabilidades.add(probs);
            }
        }
        return Probabilidades;
    }

    public static int[][] GetBestGraph(FileReader file, FileReader file2, int N) throws FileNotFoundException {
        int[][] result = new int[N][N];
        List<String[]> Grafos = new ArrayList<String[]>();
        List<String> ValoresS = new ArrayList<String>();
        try {
            BufferedReader br = new BufferedReader(file);
            BufferedReader br2 = new BufferedReader(file2);
            String str;
            String val;
            while ((str = br.readLine()) != null && (val = br2.readLine()) != null) {
                String[] lineItems = str.split("");
                String valItems = val;
                Grafos.add(lineItems);
                ValoresS.add(valItems);
            }
            br.close();
            br2.close();
        } catch (IOException e) {
            System.out.println("File not found");
        }

        List<Double> ValoresD = new ArrayList<Double>();
        for (String s : ValoresS) {
            ValoresD.add(Double.parseDouble(s));
        }

        int Menor = 2;
        /*for (int p = 0; p < ValoresD.size(); p++) {
            if (ValoresD.get(p) < ValoresD.get(Menor)) {
                Menor = p;
            }
        }*/

        int count = 0;
        for (int i1 = 0; i1 < N; i1++) {
            for (int j1 = 0; j1 < N; j1++) {
                result[i1][j1] = Integer.parseInt(Grafos.get(Menor)[count]);
                count++;
            }
        }
        count = 0;

        return result;
    }

    public static int[][] GetDataTrain(String pais1) throws Exception {
        ConverterUtils.DataSource source = new ConverterUtils.DataSource(".\\Dataset\\" + pais1 + ".arff");
        Instances newData = source.getDataSet();
        int N = newData.numAttributes();
        int trainSize = (int) Math.round(newData.numInstances());
        String[] varnames = {"PuntosFIFA", "MediaGoles", "ValorERival", "Entidad", "ValorExperto", "Resultado"};
        String[][] varvals = {{"Bajo", "Medio", "Alto"},
        {"Bajo", "Alto"},
        {"Bajo", "Alto"},
        {"UEFA", "CONMEBOL", "Otros"},
        {"Bajo", "Medio", "Alto"},
        {"Derrota", "Empate", "Victoria"}};

        Instances train = new Instances(newData, 0, trainSize);
        int[][] datatrain = new int[train.numInstances()][train.numAttributes()];
        for (int i = 0; i < train.numInstances(); i++) {
            for (int j = 0; j < train.numAttributes(); j++) {
                datatrain[i][j] = (int) train.instance(i).value(j);
            }
        }
        return datatrain;
    }

    public static int[] Getcard(String pais1) throws Exception {
        ConverterUtils.DataSource source = new ConverterUtils.DataSource(".\\Dataset\\" + pais1 + ".arff");
        Instances newData = source.getDataSet();
        int N = newData.numAttributes();
        int[] card = new int[newData.numAttributes()];
        for (int i = 0; i < newData.numAttributes(); i++) {
            card[i] = newData.attribute(i).numValues();
        }
        return card;
    }

    public static double[] GetInferencia(String pais2, int[][] result, int[] card, List<double[]> Probabilidades, int N) throws Exception {

        ConverterUtils.DataSource source = new ConverterUtils.DataSource(".\\Contrincante\\" + pais2 + ".arff");
        Instances newData2 = source.getDataSet();
        int[][] data2 = new int[N][N];
        int N3 = newData2.numAttributes();
        int[] valorTest = new int[N3];

        for (int i2 = 0; i2 < newData2.numInstances(); i2++) {
            for (int j2 = 0; j2 < newData2.numAttributes(); j2++) {
                data2[i2][j2] = (int) newData2.instance(i2).value(j2);
            }
        }

        int indice = 0;
        List<Double> inferr = new ArrayList<Double>();

        for (int ief = 0; ief < newData2.numInstances(); ief++) {

            List<Double> prv = new ArrayList<Double>();
            int vall = 0;
            for (int ic = 0; ic < newData2.numAttributes(); ic++) {
                valorTest[ic] = data2[ief][ic];
            }
            for (int vt = 0; vt < valorTest.length; vt++) {
                List<Integer> padres = getPadres(vt, result);
                int[] vars = union(vt, padres);
                int indicex = 0;

                if (vars.length == 1) {
                    indicex = valorTest[vars[0]];
                } else if (vars.length == 2) {
                    indicex = card[vars[0]] * valorTest[vars[1]]
                            + valorTest[vars[0]];
                } else if (vars.length == 3) {
                    indicex = card[vars[0]] * card[vars[1]] * valorTest[vars[2]]
                            + card[vars[0]] * valorTest[vars[1]]
                            + valorTest[vars[0]];
                } else if (vars.length == 4) {
                    indicex = card[vars[0]] * card[vars[1]] * card[vars[2]] * valorTest[vars[3]]
                            + card[vars[0]] * card[vars[1]] * valorTest[vars[2]]
                            + card[vars[0]] * valorTest[vars[1]]
                            + valorTest[vars[0]];
                } else if (vars.length == 5) {
                    indicex = card[vars[0]] * card[vars[1]] * card[vars[2]] * card[vars[3]] * valorTest[vars[4]]
                            + card[vars[0]] * card[vars[1]] * card[vars[2]] * valorTest[vars[3]]
                            + card[vars[0]] * card[vars[1]] * valorTest[vars[2]]
                            + card[vars[0]] * valorTest[vars[1]]
                            + valorTest[vars[0]];
                } else if (vars.length == 6) {
                    indicex = card[vars[0]] * card[vars[1]] * card[vars[2]] * card[vars[3]] * card[vars[4]] * valorTest[vars[5]]
                            + card[vars[0]] * card[vars[1]] * card[vars[2]] * card[vars[3]] * valorTest[vars[4]]
                            + card[vars[0]] * card[vars[1]] * card[vars[2]] * valorTest[vars[3]]
                            + card[vars[0]] * card[vars[1]] * valorTest[vars[2]]
                            + card[vars[0]] * valorTest[vars[1]]
                            + valorTest[vars[0]];
                }

                prv.add(Probabilidades.get(vt)[indicex]);
            }

            double total = 1;
            for (int iff = 0; iff < prv.size(); iff++) {
                total *= prv.get(iff);
            }
            inferr.add(total);
        }
        double sum = 0;
        for (int ii = 0; ii < inferr.size(); ii++) {
            sum += inferr.get(ii);
        }
        double[] valuee = new double[3];
        for (int ii = 0; ii < inferr.size(); ii++) {
            valuee[ii] = inferr.get(ii) / sum;
        }
        return valuee;
    }

    public static double[] Prediccion(String pais1, String pais2) throws Exception {
        double[] Prediccion = new double[3];
        int N = 6;
        int[][] datatrain = GetDataTrain(pais1);
        int[] card = Getcard(pais1);
        int[][] result = new int[N][N];
        int[][] data2 = new int[N][N];
        int count = 0;

        List<double[]> Probabilidades = new ArrayList<>();
        FileReader file = new FileReader(".\\GrafosMedida\\" + pais1 + "Grafo.txt");
        FileReader file2 = new FileReader(".\\GrafosMedida\\" + pais1 + "Medida.txt");

        result = GetBestGraph(file, file2, N);
        Probabilidades = GetProbabilidades(result, datatrain, N, card);
        Prediccion = GetInferencia(pais2, result, card, Probabilidades, N);

        return Prediccion;
    }

    public static void Resultado(String pais1, String pais2) throws Exception {

        double[] aVb = Prediccion(pais1, pais2);
        double[] bVa = Prediccion(pais2, pais1);
        String[] GEP = {"Victoria", "Empate", "Derrota"};

        
        System.out.println("************");
        System.out.println(pais1 + " VS " + pais2);
        for (int f = 0; f < aVb.length; f++) {
            System.out.println(aVb[f]);
        }
        System.out.println("************");

        System.out.println(pais2 + " VS " + pais1);
        for (int f2 = 0; f2 < aVb.length; f2++) {
            System.out.println(bVa[f2]);
        }

        int indexMayor = 0;
        for (int p = 0; p < aVb.length; p++) {
            if (aVb[p] > aVb[indexMayor]) {
                indexMayor = p;
            }
        }

        int indexMayor2 = 0;
        for (int p1 = 0; p1 < bVa.length; p1++) {
            if (bVa[p1] > bVa[indexMayor2]) {
                indexMayor2 = p1;
            }
        }

        System.out.println("************");

        if (aVb[indexMayor] > bVa[indexMayor2]) {
            System.out.println(GEP[indexMayor] + " de " + pais1 + " ante " + pais2 + " en un " + aVb[indexMayor] * 100 + "%");
        } else {
            System.out.println(GEP[indexMayor2] + " de " + pais2 + " ante " + pais1 + " en un " + bVa[indexMayor2] * 100 + "%");
        }
    }

    public static void CalculoGrafoValor() throws FileNotFoundException, Exception {

        String archivos = System.getProperty("user.dir") + "/DataSet";

        String carpetaM = System.getProperty("user.dir") + "\\Medidas\\";

        String carpetaG = System.getProperty("user.dir") + "\\Grafos\\";

        final File folder = new File(archivos);

        for (final File fileEntry : folder.listFiles()) {
            String archivo = fileEntry.getAbsolutePath();
            String nombre = fileEntry.getName().replace(".arff", "");
            ConverterUtils.DataSource source = new ConverterUtils.DataSource(archivo);
            Instances newData = source.getDataSet();
            PrintStream ps1 = new PrintStream(new File(carpetaM + nombre + "Medida.txt"));
            PrintStream ps2 = new PrintStream(new File(carpetaG + nombre + "Grafo.txt"));
            int N = newData.numAttributes();

            int trainSize = (int) Math.round(newData.numInstances());
            int testSize = newData.numInstances() - trainSize;

            Instances train = new Instances(newData, 0, trainSize);
            Instances test = new Instances(newData, trainSize, testSize);

            String[] varnames = {"PuntosFIFA", "MediaGoles", "ValorERival", "Entidad", "ValorExperto", "Resultado"};
            String[][] varvals = {{"Bajo", "Medio", "Alto"},
            {"Bajo", "Alto"},
            {"Bajo", "Alto"},
            {"UEFA", "CONMEBOL", "Otros"},
            {"Bajo", "Medio", "Alto"},
            {"Derrota", "Empate", "Victoria"}};

            int[][] data = new int[newData.numInstances()][newData.numAttributes()];
            for (int i = 0; i < newData.numInstances(); i++) {
                for (int j = 0; j < newData.numAttributes(); j++) {
                    data[i][j] = (int) newData.instance(i).value(j);
                }
            }

            int[][] datatrain = new int[train.numInstances()][train.numAttributes()];
            for (int i = 0; i < train.numInstances(); i++) {
                for (int j = 0; j < train.numAttributes(); j++) {
                    datatrain[i][j] = (int) train.instance(i).value(j);
                }
            }

            int[][] datatest = new int[test.numInstances()][test.numAttributes()];
            for (int i = 0; i < test.numInstances(); i++) {
                for (int j = 0; j < test.numAttributes(); j++) {
                    datatest[i][j] = (int) test.instance(i).value(j);
                }
            }

            int[] card = new int[newData.numAttributes()];
            for (int i = 0; i < newData.numAttributes(); i++) {
                card[i] = newData.attribute(i).numValues();
            }

            int k, i, j, w, p, m, q, z, x, c;
            long n = (int) Math.pow(2, ((Math.pow(N, 2)) - N));
            int nAristas = (int) (Math.pow(N, 2)) - N;
            int[] temp;

            int count = 0;

            ////////////////////////////////////////////// FUERZA BRUTA //////////////////////////////////////////////////
            int[][] resultFinalEntropiaFB = new int[N][N];
            int[][] resultTestEntorpiaFB = new int[N][N];
            int[][] resultFinalAkaikeFB = new int[N][N];
            int[][] resultTestAkaikeFB = new int[N][N];
            int[][] resultTestMDLFB = new int[N][N];
            int[][] resultFinalMDLFB = new int[N][N];
            int kFB = 0;
            double eTestFB;
            double eFinalFB = 0.0;
            double akaikeTestFB;
            double akaikeFinalFB = 0.0;
            double mdlTestFB;
            double mdlFinalFB = 0.0;

            Scanner s = new Scanner(System.in);
            File f = new File("Grafos.txt");
            Scanner numScan = new Scanner(f);
            String line;
            while (numScan.hasNext()) {
                line = numScan.nextLine();
                List<String[]> lineas = new ArrayList<String[]>();

                String[] lineItems = line.split(",");
                lineas.add(lineItems);
                for (k = 0; k < lineas.size(); k++) {
                    eTestFB = 0.0;
                    akaikeTestFB = 0.0;
                    mdlTestFB = 0.0;
                    for (int i1 = 0; i1 < N; i1++) {
                        for (int j1 = 0; j1 < N; j1++) {
                            resultTestEntorpiaFB[i1][j1] = Integer.parseInt(lineas.get(k)[count]);
                            resultTestAkaikeFB[i1][j1] = Integer.parseInt(lineas.get(k)[count]);
                            resultTestMDLFB[i1][j1] = Integer.parseInt(lineas.get(k)[count]);
                            count++;
                        }
                    }
                    count = 0;
                    eTestFB += GetEntropia(resultTestEntorpiaFB, card, datatrain, N);

                    if (eTestFB < eFinalFB) {
                        eFinalFB = eTestFB;
                        for (int i1 = 0; i1 < N; i1++) {
                            for (int j1 = 0; j1 < N; j1++) {
                                resultFinalEntropiaFB[i1][j1] = resultTestEntorpiaFB[i1][j1];
                            }
                        }
                    }

                    kFB = GetAkaike(GetVars(N, resultTestAkaikeFB), card);
                    akaikeTestFB = eTestFB - kFB;

                    if (akaikeTestFB < akaikeFinalFB) {
                        akaikeFinalFB = akaikeTestFB;
                        for (int i1 = 0; i1 < N; i1++) {
                            for (int j1 = 0; j1 < N; j1++) {
                                resultFinalAkaikeFB[i1][j1] = resultTestAkaikeFB[i1][j1];
                            }
                        }
                    }

                    mdlTestFB = eTestFB - ((kFB / 2) * Math.log10(train.numInstances()));

                    if (mdlTestFB < mdlFinalFB) {
                        mdlFinalFB = mdlTestFB;
                        for (int i1 = 0; i1 < N; i1++) {
                            for (int j1 = 0; j1 < N; j1++) {
                                resultFinalMDLFB[i1][j1] = resultTestMDLFB[i1][j1];
                            }
                        }
                    }
                }
            }
            ////////////////////////////////////////////// K2 //////////////////////////////////////////////////
            int[][] resultFinalEntropiaK2 = new int[N][N];
            int[][] resultTestEntorpiaK2 = new int[N][N];

            int[][] resultFinalAkaikeK2 = new int[N][N];
            int[][] resultTestAkaikeK2 = new int[N][N];
            int[][] o = new int[N][N];
            int[][] l = new int[N][N];

            int[][] resultFinalMDLFK2 = new int[N][N];
            int[][] resultTestMDLK2 = new int[N][N];

            double eTest1k2;
            double eTest2k2;
            double eFinalk2 = 0.0;
            double akaikeTest1k2;
            double akaikeTest2k2;
            double akaikeFinalk2 = 0.0;
            double mdlTest1k2;
            double mdlTest2k2;
            double mdlFinalk2 = 0.0;
            double k1 = 0.0;
            double k2 = 0.0;
            double a;
            double b;
            double a1;
            double b1;
            System.out.println("Analisando K2");
            //Entropia
            for (int g = 0; g < resultFinalEntropiaK2.length; g++) {
                for (int d = 0; d < resultFinalEntropiaK2.length; d++) {
                    eTest1k2 = 0.0;
                    eTest2k2 = 0.0;
                    if (d > g) {
                        resultTestEntorpiaK2[g][d] = 1;

                        eTest1k2 += GetEntropia(resultFinalEntropiaK2, card, datatrain, N);
                        eTest2k2 += GetEntropia(resultTestEntorpiaK2, card, datatrain, N);

                        if (eTest1k2 > eTest2k2) {
                            for (int i1 = 0; i1 < N; i1++) {
                                for (int j1 = 0; j1 < N; j1++) {
                                    resultFinalEntropiaK2[i1][j1] = resultTestEntorpiaK2[i1][j1];
                                }
                            }
                        } else {
                            resultTestEntorpiaK2[g][d] = 0;
                        }

                        eFinalk2 = eTest1k2;
                    }
                }
            }
            //Akaike
            for (int g = 0; g < resultFinalAkaikeK2.length; g++) {
                for (int d = 0; d < resultFinalAkaikeK2.length; d++) {

                    a = 0.0;
                    b = 0.0;
                    akaikeTest1k2 = 0.0;
                    akaikeTest2k2 = 0.0;
                    if (d > g) {
                        resultTestAkaikeK2[g][d] = 1;
                        a += GetEntropia(resultFinalAkaikeK2, card, datatrain, N);
                        b += GetEntropia(resultTestAkaikeK2, card, datatrain, N);

                        k1 = GetAkaike(GetVars(N, resultFinalAkaikeK2), card);
                        k2 = GetAkaike(GetVars(N, resultTestAkaikeK2), card);

                        akaikeTest1k2 = a - k1;
                        akaikeTest2k2 = b - k2;

                        if (akaikeTest1k2 > akaikeTest2k2) {

                            for (int i1 = 0; i1 < N; i1++) {
                                for (int j1 = 0; j1 < N; j1++) {
                                    resultFinalAkaikeK2[i1][j1] = resultTestAkaikeK2[i1][j1];
                                }
                            }

                        } else {
                            resultTestAkaikeK2[g][d] = 0;
                        }

                        akaikeFinalk2 = akaikeTest2k2;
                    }
                }
            }
            //MDL
            for (int g = 0; g < resultFinalMDLFK2.length; g++) {
                for (int d = 0; d < resultFinalMDLFK2.length; d++) {

                    a1 = 0.0;
                    b1 = 0.0;
                    mdlTest1k2 = 0.0;
                    mdlTest2k2 = 0.0;
                    if (d > g) {
                        resultTestMDLK2[g][d] = 1;
                        a1 += GetEntropia(resultFinalMDLFK2, card, datatrain, N);
                        b1 += GetEntropia(resultTestMDLK2, card, datatrain, N);

                        k1 = GetAkaike(GetVars(N, resultFinalMDLFK2), card);
                        k2 = GetAkaike(GetVars(N, resultTestMDLK2), card);

                        mdlTest1k2 = a1 - ((k1 / 2) * Math.log10(train.numInstances()));
                        mdlTest2k2 = b1 - ((k2 / 2) * Math.log10(train.numInstances()));

                        if (mdlTest1k2 > mdlTest2k2) {
                            for (int i1 = 0; i1 < N; i1++) {
                                for (int j1 = 0; j1 < N; j1++) {
                                    resultFinalMDLFK2[i1][j1] = resultTestMDLK2[i1][j1];
                                }
                            }
                        } else {
                            resultTestMDLK2[g][d] = 0;
                        }
                        mdlFinalk2 = mdlTest2k2;
                    }
                }
            }

            ////////////////////////////////////////////// HC //////////////////////////////////////////////////  
            temp = HipotesisMatrix(1, nAristas);
            //Entropia
            int[][] resultEntropyHC = new int[N][N];
            k = 0;
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    if (i == j) {
                        resultEntropyHC[i][j] = 0;
                    } else {
                        resultEntropyHC[i][j] = temp[k];
                        k++;
                    }
                }
            }
            Optional<Double> MaxEntopyHC = null;

            int[][] result2E = new int[N][N];
            int[][] result3E = new int[N][N];
            int[][] result4E = new int[N][N];

            TF_Topicos dc2E = new TF_Topicos(result2E);
            TF_Topicos dc3E = new TF_Topicos(result3E);
            TF_Topicos dc4E = new TF_Topicos(result4E);

            int numb = 0;

            while (numb < 500) {

                List<int[][]> GrafosHC = new ArrayList<int[][]>();
                GrafosHC = GetHC(resultEntropyHC, result2E, result3E, result4E, dc4E, dc2E, dc3E);

                List<Double> Entropias = new ArrayList<Double>();

                for (int hc2 = 0; hc2 < GrafosHC.size(); hc2++) {
                    double e = 0.0;
                    e += GetEntropia(GrafosHC.get(hc2), card, datatrain, N);
                    Entropias.add(e);

                }

                MaxEntopyHC = Entropias.stream().reduce(Double::min);
                resultEntropyHC = GrafosHC.get(Entropias.indexOf(MaxEntopyHC.get()));
                numb++;
            }
            //Akaike
            int[][] resultAkaikeHC = new int[N][N];
            k = 0;
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    if (i == j) {
                        resultAkaikeHC[i][j] = 0;
                    } else {
                        resultAkaikeHC[i][j] = temp[k];
                        k++;
                    }
                }
            }
            Optional<Double> MaxAkakiHC = null;

            int[][] result2A = new int[N][N];
            int[][] result3A = new int[N][N];
            int[][] result4A = new int[N][N];

            TF_Topicos dc2A = new TF_Topicos(result2A);
            TF_Topicos dc3A = new TF_Topicos(result3A);
            TF_Topicos dc4A = new TF_Topicos(result4A);
            int numb2 = 0;
            while (numb2 < 500) {
                List<int[][]> GrafosHCaka = new ArrayList<int[][]>();
                GrafosHCaka = GetHC(resultAkaikeHC, result2E, result3E, result4E, dc4E, dc2E, dc3E);
                List<Double> Akaike = new ArrayList<Double>();

                for (int hc2 = 0; hc2 < GrafosHCaka.size(); hc2++) {
                    double e1 = 0.0;
                    double ka = 0.0;
                    double akaikeHC = 0.0;
                    e1 += GetEntropia(GrafosHCaka.get(hc2), card, datatrain, N);

                    ka = GetAkaike(GetVars(N, GrafosHCaka.get(hc2)), card);

                    akaikeHC = e1 - ka;
                    Akaike.add(akaikeHC);
                }
                MaxAkakiHC = Akaike.stream().reduce(Double::min);
                resultAkaikeHC = GrafosHCaka.get(Akaike.indexOf(MaxAkakiHC.get()));
                numb2++;
            }
            //MDL
            int[][] resultMdlHC = new int[N][N];
            k = 0;
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    if (i == j) {
                        resultMdlHC[i][j] = 0;
                    } else {
                        resultMdlHC[i][j] = temp[k];
                        k++;
                    }
                }
            }
            Optional<Double> MaxMdlHC = null;
            int[][] result2M = new int[N][N];
            int[][] result3M = new int[N][N];
            int[][] result4M = new int[N][N];

            TF_Topicos dc2M = new TF_Topicos(result2M);
            TF_Topicos dc3M = new TF_Topicos(result3M);
            TF_Topicos dc4M = new TF_Topicos(result4M);
            int numb3 = 0;

            while (numb3 < 500) {
                List<int[][]> GrafosHCmdl = new ArrayList<int[][]>();
                GrafosHCmdl = GetHC(resultMdlHC, result2E, result3E, result4E, dc4E, dc2E, dc3E);
                List<Double> Mdl = new ArrayList<Double>();

                for (int hc2 = 0; hc2 < GrafosHCmdl.size(); hc2++) {
                    double e3 = 0.0;
                    double ka = 0.0;
                    double mdl = 0.0;
                    e3 += GetEntropia(GrafosHCmdl.get(hc2), card, datatrain, N);
                    ka = GetAkaike(GetVars(N, GrafosHCmdl.get(hc2)), card);

                    mdl = e3 - ((ka / 2) * Math.log10(train.numInstances()));
                    Mdl.add(mdl);
                }
                MaxMdlHC = Mdl.stream().reduce(Double::min);
                resultMdlHC = GrafosHCmdl.get(Mdl.indexOf(MaxMdlHC.get()));
                numb3++;

            }
            ps1.println(eFinalFB);

            for (int i8 = 0; i8 < resultFinalEntropiaFB.length; i8++) {
                for (int j5 = 0; j5 < resultFinalEntropiaFB.length; j5++) {
                    ps2.print(resultFinalEntropiaFB[i8][j5]);
                }
            }
            ps2.println();

            ps1.println(akaikeFinalFB);

            for (int i8 = 0; i8 < resultFinalAkaikeFB.length; i8++) {
                for (int j5 = 0; j5 < resultFinalAkaikeFB.length; j5++) {
                    ps2.print(resultFinalAkaikeFB[i8][j5]);
                }
            }
            ps2.println();

            ps1.println(mdlFinalFB);

            for (int i8 = 0; i8 < resultFinalMDLFB.length; i8++) {
                for (int j5 = 0; j5 < resultFinalMDLFB.length; j5++) {
                    ps2.print(resultFinalMDLFB[i8][j5]);
                }
            }
            ps2.println();

            ps1.println(eFinalk2);

            for (int i8 = 0; i8 < resultFinalEntropiaK2.length; i8++) {
                for (int j5 = 0; j5 < resultFinalEntropiaK2.length; j5++) {
                    ps2.print(resultFinalEntropiaK2[i8][j5]);
                }
            }
            ps2.println();

            ps1.println(akaikeFinalk2);

            for (int i8 = 0; i8 < resultFinalAkaikeK2.length; i8++) {
                for (int j5 = 0; j5 < resultFinalAkaikeK2.length; j5++) {
                    ps2.print(resultFinalAkaikeK2[i8][j5]);
                }
            }
            ps2.println();

            ps1.println(mdlFinalk2);

            for (int i8 = 0; i8 < resultFinalMDLFK2.length; i8++) {
                for (int j5 = 0; j5 < resultFinalMDLFK2.length; j5++) {
                    ps2.print(resultFinalMDLFK2[i8][j5]);
                }
            }

            ps1.println(MaxEntopyHC.get());

            for (int i8 = 0; i8 < resultEntropyHC.length; i8++) {
                for (int j5 = 0; j5 < resultEntropyHC.length; j5++) {
                    ps2.print(resultEntropyHC[i8][j5]);
                }
            }
            ps2.println();

            ps1.println(MaxAkakiHC.get());

            for (int i8 = 0; i8 < resultAkaikeHC.length; i8++) {
                for (int j5 = 0; j5 < resultAkaikeHC.length; j5++) {
                    ps2.print(resultAkaikeHC[i8][j5]);
                }
            }
            ps2.println();

            ps1.println(MaxMdlHC.get());

            for (int i8 = 0; i8 < resultMdlHC.length; i8++) {
                for (int j5 = 0; j5 < resultMdlHC.length; j5++) {
                    ps2.print(resultMdlHC[i8][j5]);
                }
            }
        }
    }

    public static void main(String[] args) throws Exception {

        Scanner S = new Scanner(System.in);
        String pais1;
        String pais2;
        pais1 = S.next();
        pais2 = S.next();

        Resultado(pais1, pais2);

    }
}
