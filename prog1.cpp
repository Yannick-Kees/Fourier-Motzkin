// Programmierabgabe 1
// George Tyriard & Yannick Kees
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>

using namespace std;

struct Matrix_Vector_Set // For the output from Fourier Motzkin and reconstruction
{   vector<vector<double>> Matrix;
    vector<double> Vector;
    set<int, greater<int> > Set;
};

/*
-------------------------------------------------------------------
        Some vector stuff
-------------------------------------------------------------------
*/
vector<double> get_row(vector<vector<double>> A, int j){
    // returns the j-th row of A as vector
    vector<double> sol(A[0].size());
    for (int i=0; i<A[0].size(); i++){
        sol[i] = A[j][i];
    }
    return sol;
}
double dot_product(vector<double> M, vector<double> x){
    // returns dot product of M and x
    double sum = 0.0;
    for (int i=0; i<x.size(); i++){
        sum += M[i]*x[i];
    }
    return sum;
}
void vector_print(vector<double> x){
    // print out vector
    for (int i = 0; i<x.size(); i++){
        cout << x[i] << "\t";
    }
    cout << "\n";
    return;
}
double get_min(vector<double> v  ){
    // returns minimum value of a vector
    double min = v[0];
    for (int i=0; i<v.size(); i++){
        if (v[i]<min) min = v[i];
    }
    return min;
}
double get_max(vector<double> v  ){
    // returns maximum value of a vector
    double max = v[0];
    for (int i=0; i<v.size(); i++){
        if (v[i]>max) max = v[i];
    }
    return max;
}
/*
-------------------------------------------------------------------
        Fourier Motzkin
-------------------------------------------------------------------
*/
vector<vector<int>> construct_p(set<int, greater<int> > Z, set<int, greater<int> > N, set<int, greater<int> > P){
        // returns an enumeration of Z \cup (N \times P)
        int r = Z.size() + ( N.size() * P.size());
        vector<vector<int>> p(r, vector<int>(2));

        int i = 0;

        // set p = Z \cup (N x P)
        for (set<int, greater<int> >::iterator itr_z = Z.begin(); itr_z != Z.end(); itr_z++)
        {   // add Z
            p[i][0] = *itr_z ;
            p[i][1] = -1; // if an element is in Z, set the second coordinate -1
            i++;
        }
        for (set<int, greater<int> >::iterator itr_n = N.begin(); itr_n != N.end(); itr_n++)
        {
            for (set<int, greater<int> >::iterator itr_p = P.begin(); itr_p != P.end(); itr_p++){
                // add N times P
                p[i][0] = *itr_n ;
                p[i][1] = *itr_p;
                i++;
            }
        }
        return p;
        
    }


Matrix_Vector_Set Fourier_motzkin(std::vector<std::vector<double>> A, std::vector<double> b, int j){
    // Eliminate the j-th variable of Ax <= b
    int n = A[0].size();

    set<int, greater<int> > Z; // zero indices
    set<int, greater<int> > N; // negative indices
    set<int, greater<int> > P; // positive indices
    
    for (int i =0; i< A.size(); i++  )
    {
        if (A[i][j] > 0) {
            P.insert(i); 
        }
        if (A[i][j] < 0) {
            N.insert(i);
        }
        if (A[i][j] == 0) {
            Z.insert(i);
        }
    }
    vector<vector<int>> p = construct_p(Z, N, P);
    int r = p.size();
    vector<double> d(r);
    vector<vector<double>> D(r, vector<double>(n  ));

    for (int i=0; i<r; i++){
        if (p[i][1]== -1){ // if p[i] in Z
            // D_i. <-A_p(i).
            for(int k=0; k<n; k++){
                D[i][k] = A[p[i][0]][k];
            }
            d[i] = b[p[i][0]];

        } else{ // if p(i) in N x P
            int s = p[i][0];
            int t = p[i][1];
            // D_i <- a_tj * A_s - a_sj * A_t
            for(int k=0; k<n; k++){
                D[i][k] = A[t][j] * A[s][k] - A[s][j] * A[t][k];
            }
            d[i] = A[t][j] * b[s] - A[s][j] * b[t];
        }
    }
    Matrix_Vector_Set output; // we cannot return 2 outputs in C++, so use the struct
    output.Matrix = D;
    output.Vector = d;
    output.Set = P;
    return output;
}

/*
-------------------------------------------------------------------
        Get feasable vector
-------------------------------------------------------------------
*/

bool is_valid(vector<double> d){
    // by Fourier-Motzkin there exists a solution, if in the last iteration, all values of the rhs are positive
    for(int i=0; i<d.size(); i++){
        if ( d[i] < 0) return false;
    }
    return true;
}


void calculate_feasable_vector(vector<Matrix_Vector_Set> V){
    vector<double> x(V[0].Matrix[0].size(), 0.0);
    for(int i = V.size()-1; i>= 0; i--){
        if(V[i].Set.size()>0){
            // If U is not empty we choose the minimum over these positiv inequalitys
            vector<double> b_i_a_ix(V[i].Set.size());
            int k = 0;

            for (set<int, greater<int> >::iterator itr = V[i].Set.begin(); itr != V[i].Set.end(); itr++){
                b_i_a_ix[k] = ( 1.0/V[i].Matrix[*itr][i] ) *  (V[i].Vector[*itr] - dot_product(get_row(V[i].Matrix, *itr), x));
                k++;
            }
            
            x[i] = get_min(b_i_a_ix);
        } else{// what do we do, if we only have non-positiv a_ik values
            vector<double> b_i_a_ix(V[i].Matrix.size(), -1000.0);
            for (int k=0; k<V[i].Matrix.size(); k++  ){
                if (V[i].Matrix[k][i] != 0){
                b_i_a_ix[k] = ( 1.0/V[i].Matrix[k][i] ) *  V[i].Vector[k] - dot_product(get_row(V[i].Matrix, k), x);
                k++;
                }
            }
            x[i] = get_max(b_i_a_ix);
        }   
    }
    vector_print(x);
    return;
}
/*
-------------------------------------------------------------------
        Main function
-------------------------------------------------------------------
*/
int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " FILENAME\n";
        return 1;
    }
     /*
    -------------------------------------------------------------------
            Get Input
    -------------------------------------------------------------------
    */

    ifstream input(argv[1]);

    if(!input){
        cerr << "Fehler beim Öffnen der Datei " << argv[1] << "\n";
        return 1;
    }

    int m = 0; // #rows
    int n = 0; // #columns

    input >> m >> n;

    vector<double> c(n); // get c
    for (int i = 0; i< n; i++){
        input >> c[i];
    }
    
    vector<double> b(m) ; // get b
    for (int i = 0; i< m; i++){
        input >> b[i];
    }
    vector<vector<double>> A(m, vector<double>(n)); // get A
    for (int i = 0; i< m; i++){
        for (int j = 0; j< n; j++){
            input >> A[i][j];
        }
    }

     /*
    -------------------------------------------------------------------
            Process data
    -------------------------------------------------------------------
    */
    // We store everything we need for the ouput vectors in this array
    vector<Matrix_Vector_Set> everything_important(n);

    for(int j=0; j<n; j++ ){
        // Eliminate all Variables and save the produced part-solutions
        Matrix_Vector_Set construct = Fourier_motzkin(A, b, j);
        everything_important[j].Matrix = A;
        everything_important[j].Vector = b;
        everything_important[j].Set = construct.Set;
        A = construct.Matrix;
        b = construct.Vector;
    }

    if (is_valid(b)){
        calculate_feasable_vector(everything_important);
    } else{
        cout << "empty ";
        A = everything_important[0].Matrix;
        b = everything_important[0].Vector;
        /*
        -------------------------------------------------------------------
                Get the certificate vector
        -------------------------------------------------------------------
        
        We solve
                ( A^T )         ( 0 )
                ( -A^T )  u <=  ( 0 )
                ( -id )         ( 0 )
                ( b^t )         ( -1 )

        Note that this is eqivalent to the one in the Lemma of Farkas.
        We know that if A^Ty=0 then A^t * lambda * y = 0 for every lambda in R, and since y>=0 we know that
        b^Ty < 0   <=>   b^Ty< -1 
        So we can do what we did before on this LP.
        */
        vector<vector<double>> AA(2*n+m+1, vector<double>(m));
        //initialise the new lhs 
            for(int i=0; i< A.size(); i++){
                for(int j=0; j<A[0].size(); j++){
                    AA[j][i] = A[i][j];
                    AA[j+n][i] = (-1.0) * A[i][j];
                }
            }
            
            for(int i=0; i< m; i++){
                AA[i+2*n][i] = -1;
            }
            for(int i=0; i<m; i++){
                AA[2*n+m][i] = b[i];
            }


        vector<double> bb(2*n+m+1, 0.0);
        //initialiste the new rhs
            bb[2*n+m] = -1;

        // this is now the same as before..
        vector<Matrix_Vector_Set> everything_important2(m);
        
        for(int j=0; j<m; j++ ){
            // Eliminate all Variables and save the produced part-solutions
            Matrix_Vector_Set construct2 = Fourier_motzkin(AA, bb, j);
            everything_important2[j].Matrix = AA;
            everything_important2[j].Vector = bb;
            everything_important2[j].Set = construct2.Set;
            AA = construct2.Matrix;
            bb = construct2.Vector;
        }

     
        if (is_valid(bb)){
            // the feasable vector for the new LP is our certificate vector
            calculate_feasable_vector(everything_important2);
        } else{
            printf("Fehler");
        }
        
    }
}