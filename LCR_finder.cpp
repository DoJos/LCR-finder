#include <iostream>
#include <chrono>
#include <utility>
#include <math.h>
#include <tuple>
#include <algorithm>
#include <vector>
#include <stdlib.h>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/score.h>
#include <seqan/align.h>
#include <seqan/score/score_matrix_data.h>
#include <seqan/seeds.h>
#include <seqan/find.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/file.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>

//valgrind --tool=callgrind ./my_proj ...


/*
 * TODOS
 * 
 * ES GIBT EIN SEED EXTENSION ALGORITHMUS IN SEQAN
 * 
 * REPEAT Matrix anpassen
 * 2. kondition zum vertex erstellen einbauen
 *
 * finding longest path ueberpruefen
 * extending longest-path intervals einbauen
 *
 * argumente als parameter übergeben
 *
 *
 *
 * Git fork erstellen
 */

using namespace seqan;
using namespace std::chrono;
using namespace std;

typedef unsigned int TCargo;
typedef Graph<Directed<TCargo> > TGraph;
typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
typedef Size<TGraph>::Type TSize;
typedef std::chrono::high_resolution_clock::time_point TimeVar;
typedef double TScoreValue;

//complexity cutoff value for extending longest path intervals
//the cutoff value was calculated by Li & Kahveci by randomly sampling sequences from Swissprot that contain repeat regions
// (Li & Kahveci, Bioinformatics 2006)
double comCut = 3.7491631225622157;

//alphabet size, 20 AminoAcids, String<AminoAcid> contains 27 letters (wildcards..)
int const alphSize = 20;


//
//reading the R-NR-Matrices from a file
Score<TScoreValue, ScoreMatrix<AminoAcid, Default> > readRNRMatrix ( String<char> RNRMatrixPath )
{
    Score<TScoreValue, ScoreMatrix<AminoAcid, Default> > RNRMatrix;
    loadScoreMatrix ( RNRMatrix, toCString ( getAbsolutePath ( toCString ( RNRMatrixPath ) ) ) );
    return RNRMatrix;
}


//GLOBAL MATRICES
String<char> nonRepeatMatrixPath = "../my_project-build_Kdevelop/combinedMatricesRowByRow095NonRepeat";
String<char> repeatMatrixPath = "../my_project-build_Kdevelop/combinedMatricesRowByRow095Repeat";

Score<TScoreValue, ScoreMatrix<AminoAcid, Default> > nonRepeatMatrix = readRNRMatrix ( nonRepeatMatrixPath );
Score<TScoreValue, ScoreMatrix<AminoAcid, Default> > repeatMatrix = readRNRMatrix ( repeatMatrixPath );
const Blosum62 blosum62Matrix;
double twoLetterScoringMatrix[676][676] = { {} };




/* Initiating 2-gram String for any combination of letters in the alphabet
 * 
 */
StringSet<String<char> > initTwoGramAlphabet() {
    String<char> alphabet;
    StringSet<String<char> > alphabetTwoGram;
    resize (alphabet, 26);
    
    for (unsigned i = 0; i < 26; ++i) {
        alphabet[i] = i + 65;
    }
    
    for (unsigned i = 0; i < 26; ++i) {
        //are similar letters also a k-gram? if not, change j=0 to j=1;
        for (unsigned j = 0; j < 26; ++j) {
            String<char> twoGram = alphabet[i];
            appendValue (twoGram, alphabet[j]);
            appendValue (alphabetTwoGram, twoGram);
            clear (twoGram);
        }
    }
    
    return alphabetTwoGram;
}


/* normalize the Scoring Matrix (Blosum62) by the sum of each column
 * 
 */
template <size_t rows, size_t columns>
void normalizeScoringMatrix (double (&array)[rows][columns]) {
    
    for (unsigned i = 0; i < rows; ++i) {
        for (unsigned j = 0; j < columns; ++j) {
            array[i][j] = pow (2, array[i][j]);
        }
    }
        
        
    double sum = 0;
    
    
    for ( unsigned i = 0; i < rows; ++i ) {
        sum = 0;
        for ( unsigned j = 0; j < columns; ++j) {
            sum = sum + array[i][j];
        }
        for (unsigned j = 0; j < columns; ++j) {
            array[i][j] = (array[i][j]) / sum;
        }
    }
    
    return;
}


/* normalize the Blosum62Matrix for primary complexity measure and store it in normalizedScoringMatrix
 * TODO USELESS ?!
 */

/*
void normalizeMatrixAndCalcSimilarityVec (int rows, int columns, String<int> & wholeSeqFrequencies) {
    double normalizedScoringMatrix[rows][columns];
    StringSet<String<char> > TwoGramList = initTwoGramAlphabet();

    for ( unsigned i = 0; i < rows; ++i ) {
        for ( unsigned j = 0; j < columns; ++j ) {
            char rowLet = i+65;
            char colLet = j+65;
            // why 2^score ?
            normalizedScoringMatrix[i][j] = pow ( 2, score ( blosum62Matrix, rowLet, colLet ) );
        }
    }


    // Print calculated normalized Scoring matrix:
    /*
    for ( unsigned i = 0; i < rows; ++i ) {
        for ( unsigned j = 0; j < columns; ++j ) {
            std:: cout << normalizedScoringMatrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    

    //normalize by the sum of each column
    double sum = 0;
    for ( unsigned i = 0; i < rows; ++i ) {
        sum = 0;
        for ( unsigned j = 0; j < columns; ++j ) {
            sum = sum + normalizedScoringMatrix[i][j];
        }
        for ( unsigned j = 0; j < columns; ++j ) {
            // why divide by the sum?
            normalizedScoringMatrix[i][j] = ( normalizedScoringMatrix[i][j] ) /sum;
        }
    }

    
    // Print calculated normalized Scoring matrix:

    for ( unsigned i = 0; i < rows; ++i ) {
        for ( unsigned j = 0; j < columns; ++j ) {
            std:: cout << normalizedScoringMatrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    
    
    //initialize similarity Vector
    String<double> similarityVector;
    resize(similarityVector, 123);
    for (unsigned i = 65; i<91; ++i) {
        similarityVector[i] = 0.;
    }
    
    //calculate similarityVector
    for (unsigned i = 65; i < 91; ++i) {
        for (unsigned j = 65; j < 91; ++j) {
            similarityVector[i] = similarityVector[i] + normalizedScoringMatrix[i-65][j-65]*wholeSeqFrequencies[j];
        }
    }
  
    String<double> normalizedSimilarityVector;
    resize(normalizedSimilarityVector, 123);
    double sumOfSimilarity = 0;
    for (unsigned i = 0; i < 26; ++i) {
        if (wholeSeqFrequencies[i] != 0) {
            sumOfSimilarity += similarityVector[i+65];
        }
    }
    
    //TODO sumOfSimilarity is 0 ?
    std::cout <<"sumOfSimilarity: " << sumOfSimilarity << std::endl;
    
    //hier oder im vorherigen schritt: moeglicherweise verfälschung der egbenisse durch AminosaeureGruppen wie B, J, O etc
    
    //TODO 
    for (unsigned i = 0; i < 26; ++i) {
        normalizedSimilarityVector[i+65] = similarityVector[i+65]/sumOfSimilarity;
        //std::cout << "normalizedSimilarityVector: " << normalizedSimilarityVector[i+65] << std::endl;
    }
    
    return;
}
*/

template <typename T>
void quickSort (T & arr, int left, int right) {
    int i = left;
    int j = right;
    
    int tmp;
    int pivot = arr[(left + right) / 2];
    
    while (i <= j) {
        while (arr[i] < pivot) {
            ++i;
        }
        while (arr[j] > pivot) {
            --j;
        }
        if (i <= j) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            ++i;
            --j;
        }
    }
    
    if (left < j) {
        quickSort(arr, left, j);
    }
    if (i < right) {
        quickSort(arr, i, right);
    }
}

/* Reversing the order of an array or string e.g. GMNK -> KNMG
 * 
 */
template <typename T>
void reverseOrder (T & arr) {
    typename Value<T>::Type tmp = arr[0];
    int size = length(arr);
    
    for (int i = 0, j = size-1; i < size/2; ++i, --j) {
        tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
    
    return;
}


/* Print the constructed graph to a file named graph.dot
 * 
 */
void drawGraph ( TGraph & myGraph ) {
    std::ofstream dotFile ( "graph.dot" );
    writeRecords ( dotFile, myGraph, DotDrawing() );
    dotFile.close();
    //std::cout << myGraph << '\n';
    
    return;
}


template <typename seqAnString>
void printSeqAnString (seqAnString & str) {
    for (unsigned i = 0; i < length(str); ++i) {
        std::cout << str[i] << " ";
    }
    std::cout << '\n';
    
    return;
}


/* prints the vertex list
 * 
 */
template <typename TVertexArray>
void printVertexArray ( TVertexArray const & vertexArray ) {
    for ( int i = 0; i < length ( vertexArray ); ++i ) {
        for ( int x = 0; x< 3; ++x ) {
            std::cout <<  vertexArray[i][x] << ' ';
        }
        std::cout << '\n';
    }
    
    return;
}


// StringSet<String<int> > & vertexArray, String<AminoAcid> const & seq
//template <typename TVertexArray, typename TSeq>
/* create Edges in our Graph with specific conditions
 * cutoff = 1, tTwo = 3, tThree = 5
 */
TGraph createEdges ( TGraph myGraph, String<String<int> > & vertexArray, String<AminoAcid> const & seq ) {
    int cutoff = 1;
    int tTwo = 3;
    int tThree = 5;

    if ( length ( vertexArray ) <= 1 ) {
        return myGraph;
    }

    //set predecessors
    int dummyPos = 0;
    for (int i = 1; i < length(vertexArray); ++i) {
        addEdge (myGraph, dummyPos, i);
    }
    

    for ( int startVertex = 1; startVertex < length ( vertexArray ); ++startVertex ) {
        //std::cout << "STARTVERTEX: " << vertexArray[startVertex][1] << "," << vertexArray[startVertex][2] << std::endl;
        for ( int destinationVertex = startVertex + 1; destinationVertex < length ( vertexArray ); ++destinationVertex ) {
            //std::cout << "DESTINATIONVERTEX: " << vertexArray[destinationVertex][1] << "," << vertexArray[destinationVertex][2] << std::endl;
            std::tuple<char, char> first ( seq[vertexArray[startVertex][1]-1], seq[vertexArray[destinationVertex][1]-1] );
            std::tuple<char, char> second ( seq[vertexArray[startVertex][2]-1], seq[vertexArray[destinationVertex][2]-1] );

            if ( score ( blosum62Matrix, std::get<0> ( first ), std::get<0> ( second ) ) > cutoff && score ( blosum62Matrix, std::get<1> ( first ), std::get<1> ( second ) ) > cutoff ) {
                //std::cout << seq[vertexArray[startVertex][1]] << vertexArray[startVertex][1] << seq[vertexArray[destinationVertex][1]] << vertexArray[destinationVertex][1] << std::endl;
                //std::cout << seq[vertexArray[startVertex][2]] << vertexArray[startVertex][2] << seq[vertexArray[destinationVertex][2]] << vertexArray[destinationVertex][2] << std::endl;

                // 3 conditions to be checked before creating an edge
                if ( abs ( ( vertexArray[startVertex][2] - vertexArray[startVertex][1] ) - ( vertexArray[destinationVertex][2] - vertexArray[destinationVertex][1] ) ) <= tTwo ) {
                    if ( vertexArray[startVertex][1] <= vertexArray[destinationVertex][1] && vertexArray[destinationVertex][1] <= vertexArray[startVertex][2] && vertexArray[startVertex][2] <= vertexArray[destinationVertex][2] ) {
                        if ( vertexArray[startVertex][2] == vertexArray[destinationVertex][1] ) {
                            if ( vertexArray[startVertex][2]-vertexArray[startVertex][1] <= tThree && vertexArray[destinationVertex][2]-vertexArray[destinationVertex][1] ) {
                                addEdge ( myGraph, startVertex, destinationVertex );
                            }
                        } else {
                            addEdge ( myGraph, startVertex, destinationVertex );
                        }
                    }
                }
            }
        }
    }

    drawGraph ( myGraph );
    return myGraph;
}


// String<int> & vertex, StringSet<String<int> > & vertexArray, int const counter, int & firstPositionOfTuple, int & secondPositionOfTuple
//template <typename TVertex, typename TVertexArray>
void createVertexArray ( String<int> & vertex, String<String<int> > & vertexArray, int const counter, int & firstPositionOfTuple, int & secondPositionOfTuple, String<int> & alphabetFreq ) {
    clear (vertex);
    appendValue (vertex, counter);
    appendValue (vertex, firstPositionOfTuple);
    appendValue (vertex, secondPositionOfTuple);

    resize ( vertexArray, counter+1 );
    vertexArray[counter] = vertex;

    return;
}


/* Checking if two aminoacids are in the same window by chance
 * 
 */
template <typename TScoreValue, typename TSequenceValue, typename TSpec>
bool checkSameWindowByChance ( Score<TScoreValue, ScoreMatrix<TSequenceValue, TSpec> > const & repeatMatrix, Score<TScoreValue, ScoreMatrix<TSequenceValue, TSpec> > const & nonRepeatMatrix, char const & firstLetter, char const & secondLetter, double const & alphabetFreq, int const & winLen ) {
    double scoreRepeatMatrix = score ( repeatMatrix, firstLetter, secondLetter );
    double scoreNonRepeatMatrix = score ( nonRepeatMatrix, firstLetter, secondLetter );
    double frequency = alphabetFreq/winLen;

    //std::cout << "|"<<scoreRepeatMatrix << "-" << frequency<<"|" << "<="<< "|"  << scoreNonRepeatMatrix << "-" << frequency<< "|" << std::endl;
    if (abs(scoreRepeatMatrix - (frequency)) <= abs(scoreNonRepeatMatrix - (frequency))) {
        return true;
    }
    
    return false;
}


/* initialize the alphabet vector for calculating the frequency of characters in the sliding window
 * 
 */
String<int> initFreqVector () {
    String<int> alphabetFreq;
    
    //sized 123 because the last letter 'z' in ASCII is at position 122
    resize (alphabetFreq, 123);
    for (unsigned i = 0; i < 123; i++) {
        alphabetFreq[i] = 0;
    }

    return alphabetFreq;
}


/* set alphabet frequencies to 0
 * 
 */
void clearAlphabetFreq(String<int> & alphabetFreq) {
    for (unsigned i = 65; i < 91; i++) {
        alphabetFreq[i] = 0;
    }
    
    return;
}

/* calculate the frequency of each character in the given sequence
 * 
 */
void getAminoAcidCountInSequence (String<int> & wholeSeqFrequencies, String<AminoAcid> const & seq) {
    for (unsigned i = 0; i < length(seq); ++i) {
        char letter = seq[i];
        ++wholeSeqFrequencies[letter];
    }
    
    return;
}


//String<AminoAcid> const & seq, String<int> & vertex, StringSet<String<int> > & vertexArray
/* create a vertex for each 2 aminoacids in a window which have a higher Blosum62 score than a given cutoff (here: 1)
 * 
 */
template <typename TSeq, typename TVertex, typename TVertexArray>
TGraph fastCreateVertices ( TSeq const & seq, unsigned const & winLen, unsigned const & maxIndelTandem, unsigned const & maxIndelCryptic, TVertex & vertex, TVertexArray & vertexArray ) {
    TGraph gbaGraph;
    int cutoff = 1;
    int counter = 0;
    //dummySource :
    int firstDumPos = 0;
    int secDumPos = 0;
    
    String<int> alphabetFreq = initFreqVector();

    //calc frequencyVector - TODO: there is a calc freqVector function
    for ( int i = 0; i < winLen; ++i ) {
        char currentLetter = seq[i];
        ++alphabetFreq[currentLetter];
    }
    
    createVertexArray ( vertex, vertexArray, counter, firstDumPos, secDumPos, alphabetFreq );
    addVertex(gbaGraph);
    
    ++counter;
    
    for ( unsigned i = 0; i < winLen; ++i ) {
        char firstLetter = seq[i];
        
        for ( unsigned x = i+1; x < ( winLen ); ++x ) {
            char secondLetter = seq[x];

            if ( score ( blosum62Matrix, firstLetter, secondLetter ) > cutoff ) {

                if ( checkSameWindowByChance ( repeatMatrix, nonRepeatMatrix, firstLetter, secondLetter, alphabetFreq[firstLetter], winLen ) ==true ) {
                    int firstPositionOfTuple = i+1;
                    int secondPositionOfTuple = x+1;
                    createVertexArray ( vertex, vertexArray, counter, firstPositionOfTuple, secondPositionOfTuple, alphabetFreq );
                    addVertex ( gbaGraph );

                    ++counter;
                }
            }
        }
    }

    for ( unsigned i = 1; i < ( length ( seq )-winLen ) + 1; ++i ) {
        //for calculating the frequencies of aminoacids in the same window
        char decrementLetter = seq[i-1];
        char incrementNextLetter = seq[i+winLen-1];

        --alphabetFreq[decrementLetter];
        ++alphabetFreq[incrementNextLetter];

        for ( int x = 0; x < winLen - 1; ++x ) {
            char firstLetter = seq[i+x];
            char secondLetter = seq[i+winLen-1];

            if ( score ( blosum62Matrix, firstLetter, secondLetter ) > cutoff ) {
                if ( checkSameWindowByChance ( repeatMatrix, nonRepeatMatrix, firstLetter, secondLetter, alphabetFreq[firstLetter], winLen ) ==true ) {
                    int firstPositionOfTuple = i+x+1;
                    int secondPositionOfTuple = i+winLen-1+1;

                    createVertexArray ( vertex, vertexArray, counter, firstPositionOfTuple, secondPositionOfTuple, alphabetFreq );
                    addVertex ( gbaGraph );

                    ++counter;
                }
            }
        }
    }

    return gbaGraph;
}


/* Calculate shannon entropy
 * 
 */
long double shanEntropy ( unsigned const & L, String<int> const & aaQuantity ) {
    long double K2KV = 0;
    
    for ( unsigned i = 65; i < 26+65; ++i ) {    
        if ( aaQuantity[i] != 0 ) {
            K2KV += -log2((long double)aaQuantity[i] / L ) * ( (long double)aaQuantity[i] / L );
        }
    }
    
    return K2KV;
}

/* initialize sequence-entropies with -1
 * 
 */
String<double> initSeqEnt(String<AminoAcid> const & seq) {
    int seqLen = length(seq);
    String<double> seqEnt;
    resize (seqEnt, seqLen);
    
    for (unsigned i = 0; i < seqLen; ++i) {
        seqEnt[i] = -1;
    }
        
    return seqEnt;
}

/* TODO
 * 
 */
int findLow (float const & hiCut, int const & i, int const & limit, String<long double> const & H) {
    int j;
    
    for (j = i; j >= limit; --j) {    
        if (H[j] == -1) {
            break;
        } else if (H[j] > hiCut) {
            break;
        }
    }
    
    return (j+1);
}

/* TODO
 * 
 */
int findHigh (float const & hiCut, int const & i, int const & limit, String<long double> const & H) {
    int j;
    
    for (j = i; j <= limit; ++j) {
        if (H[j] == -1) {
            break;
        } else if (H[j] > hiCut) {
            break;
        }
    }
    
    return (j-1);
}

/* calculates factorial
 * 
 */
long double factorial (long double n) {
    if (n == 0) {
        return 1;
    } else {
        return n * factorial (n-1);
    }
}

/* calculates the natural log of factorial of n
 * 
 */
long double lnFact (double n) {
    long double ans;
    ans = factorial(n);
    ans = log(ans);
    
    return ans;
}

/*Calculates shan Entropy per equation 3 of Wootton and Federhen (March 1993)
 * PERMUTATION (natural log)
 */
long double lnPerm (String<int> const & freqVector, int winLen) {
    long double ans;
    int i;
    
    ans = lnFact(winLen);
    
    for (i = 0; i < length(freqVector); ++i) {
        if (freqVector[i] != 0) {
            ans -= lnFact(freqVector[i]);
        }
    }
    
    return ans;
}


/* Calculates nat Log of colourings, the number of compositions in any complexity state
 * See equation 1 of Wootton and Federhen (March 1993)
 * 
 * */
double lnAss (String<int> const & freqVector) {
    long double ans;
    
    //TODO freqVector hat viele unnötige stellen, die alle 0 sind, da  das alphabet in ASCII erst bei 65 (?) beginnt
    String<int> sortedVector = freqVector;
    quickSort(sortedVector, 0, length(sortedVector)-1);
    reverseOrder(sortedVector);
    
    ans = lnFact(alphSize);
    if (sortedVector[0] == 0) {
        return ans;
    }
    
    int total = alphSize;
    int cl = 1;
    
    int svim = sortedVector[0];
    int svi = sortedVector[0];
    
    for (unsigned i = 0; i < length(sortedVector); ++i) {    
        svim = svi;
    
        if (i == alphSize) {
            ans -= lnFact(cl);
            break;
        } else if ((svi = sortedVector[i+1]) == sortedVector[i]) {
            cl++;
            continue;
        }
        
        else {
            total -= cl;
            ans -= lnFact(cl);
            
            if (svi == 0) {    
                ans -= lnFact(total);
                break;
            } else {
                cl = 1;
                continue;
            }
        }
    }
    
    return ans;
}


/* getProb(), ln(p0) in Wootton and Federhen
 * 
 */
long double getProb (int & len, String<int> const & freqVector) {
    long double ans, ans1, ans2 = 0.;
    long double totSeq = ((long double) len) * (log((long double) alphSize));
    ans1 = lnAss (freqVector);
    
    if (ans1 > -100000.0) {
        ans2 = lnPerm (freqVector, len);
    } else {
        //TODO
        std::cout << " INSERT ERROR IN FUNCTION GETPROB" << std::endl;
    }
    ans = ans1 + ans2 - totSeq;
    
    
    return ans;
}


/* trim() is trimming the subsequence by calculating the probability distributiion of complexity states
 * for a given window length.
 * 
 */
void trim(String<AminoAcid> const & tempSeq, int & leftEnd, int & rightEnd, int const winLen) {
    long double prob, minProb;
    int len;
    int lend, rend;
    int minLen;
    int maxTrim;
    
    lend = 0;
    rend = length(tempSeq) -1;
    minLen = 1;
    maxTrim = 100;
    
    if ((length(tempSeq)-maxTrim) > minLen) {
        minLen = length(tempSeq)-maxTrim;
    }
    
    minProb = 1.;
    
    //shortens the window by 1 every time
    for (len = length(tempSeq); len > minLen; --len) {
        int count = 0;
        int i = 0;
        int winStartPos = 0;
        int winEndPos = len-1;
        
        if (winEndPos <= 0) {
            break;
        }
        
        //shift window of length len over subsequence
        for (unsigned x = 0; x < length(tempSeq)-len+1; ++x) {
            //TODO in folgenderzeile winStartPos+count+2 oder +1 oder 0
            String<AminoAcid> window = infix (tempSeq, winStartPos+count, winEndPos+count+1);
            
            if (length(window) <= 1) {
                break;
            }
            
            String<int> freqVector = initFreqVector();
            getAminoAcidCountInSequence(freqVector, window);
            prob = getProb (len, freqVector);
            
            if (prob < minProb) {
                minProb = prob;
                lend = i;
                rend = len + i - 1;
            }
            
            ++count;
            ++i;
        }
    }
    leftEnd = leftEnd + lend;
    rightEnd = rightEnd - (length(tempSeq) - rend - 1);
    
    return ;
    
}


/* getSourceVertex() finds all vertices that have no incoming edges
 * 
 */
String<int> getSourceVertex ( TGraph const & gbaGraph )
{
    String<int> sourceVertexList;
    int vertexCount = numVertices ( gbaGraph );

    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it ( gbaGraph );
    for ( unsigned i = 0; i < vertexCount; ++i ) {
        if ( inDegree ( gbaGraph, i ) == 0 ) {
            appendValue ( sourceVertexList, i );
        }
    }

    return sourceVertexList;
}


/* getPaths() finds all vertex IDs of a path longer than a cutoff (3)
 * 
 */
template <typename TSpec, typename TPredecessorMap, typename TVertexDescriptor1, typename TVertexDescriptor2>
void getPaths ( Graph<TSpec> const& g,
                TPredecessorMap const& predecessor,
                TVertexDescriptor1 const source,
                TVertexDescriptor2 const v,
                String<int> & longPaths
              )
{

    if ( source == v ) {
        //std::cout << source;
        appendValue ( longPaths, source );
    } else if ( getProperty ( predecessor, v ) == getNil<typename VertexDescriptor<Graph<Directed<TCargo> > >::Type>() ) {
        //std::cout << "No path from " << source << " to " << v << " exists.";

    } else {
        getPaths ( g,predecessor, source, getProperty ( predecessor, v ), longPaths );
        //std::cout << "," << v;
        appendValue ( longPaths, v );
    }

}


/* getLongestPaths() returns the longest path in a directed graph (reversed dijkstra algorithm, edgeValue * -1)
 * 
 */
String<int> getLongestPaths ( TGraph const & gbaGraph )
{
    int edgeCount = numEdges ( gbaGraph );
    String<int> sourceVertices = getSourceVertex ( gbaGraph );
    int source = 0;
    String<int> weights;
    resize ( weights, edgeCount );
    String <int> weightMap;
    String<int> predMap;
    String<int> distMap;
    
    for ( int i = 0; i < numEdges ( gbaGraph ); ++i ) {
        weights[i] = -1;
    }

    assignEdgeMap ( weightMap, gbaGraph, weights );    
    dagShortestPath ( predMap, distMap, gbaGraph, source,weightMap );

    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it ( gbaGraph );
    String<int> longPaths;

    while ( !atEnd ( it ) ) {
        if ( getProperty ( distMap, getValue ( it ) ) <= -3 ) {
            getPaths ( gbaGraph, predMap, ( TVertexDescriptor ) source, getValue ( it ), longPaths );
            //std::cout << " (Distance: " << getProperty ( distMap, getValue ( it ) ) << ")\n";
        }
        goNext ( it );
    }

    return longPaths;
}

// ---------------------------EXECUTION TIME FUNCTION

#define duration(a) std::chrono::duration_cast<std::chrono::microseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()
template<typename F, typename... Args>
double funcTime ( F func, Args&&... args )
{
    TimeVar t1=timeNow();
    func ( std::forward<Args> ( args )... );
    return duration ( timeNow()-t1 );
}


/* maskSequence() mask the sequence at the given region by changing aminoacids to character X
 * 
 */
void maskSequence ( String<AminoAcid> & seq, int const startMaskPosition, int const endMaskPosition )
{
    String<AminoAcid> maskedSeq;
    for ( int i = startMaskPosition; i <= endMaskPosition; ++i ) {
        seq[i-1] = 'X';
    }
    
    return;
}


/* showScoringMatrix()
 * 
 */
template <typename TScoreValue, typename TSequenceValue, typename TSpec>
void showScoringMatrix ( Score<TScoreValue, ScoreMatrix<TSequenceValue, TSpec> > const & scoringScheme )
{


    for ( unsigned i=0; i<ValueSize<TSequenceValue>::VALUE; ++i ) {
        std::cout << "\t" << TSequenceValue ( i );
    }

    std::cout << std::endl;

    for ( unsigned i = 0; i < ValueSize<TSequenceValue>::VALUE; ++i ) {
        for ( unsigned j=0; j < ValueSize<TSequenceValue>::VALUE; ++j ) {
            std::cout << "\t" << score ( scoringScheme, TSequenceValue ( i ), TSequenceValue ( j ) );
        }
        std::cout << std::endl;
    }

    return;
}

struct sampLenF {
    std::string start;
    std::string end;
    std::string mean;
    std::string sd;
    
    friend std::istream& operator>>(std::istream& str, sampLenF& data) {
        std::string line;
        sampLenF tmp;
        if (std::getline(str,line)) 
        {
            std::stringstream iss(line);
            if (std::getline(iss, tmp.start, ',')       &&
                std::getline(iss, tmp.end, ' ')         &&
                std::getline(iss, tmp.mean, ':')        &&
                std::getline(iss, tmp.sd))
            {
                data.swap(tmp);
            } else {
                std::cout << "Filterprocess failed. Can't read file 'sampledLenRepPer'." << std::endl;
                //str.setstate(std::ios::failbit);
            }
        }
     
        return str;
    }
    
    void swap(sampLenF& other) {
        std::swap(start, other.start);
        std::swap(end, other.end);
        std::swap(mean, other.mean);
        std::swap(sd, other.sd);
    }
};

/* readSampleLenPercentage()
 * 
 */
String<String<double> > readSampleLenPercentage (String<char> const & pathToFile) {
    //String<char> sampleLenPerFile = getAbsolutePath (toCString(pathToFile));
    String<String<double> > sampLenPer;
    String<double> grp;
    ifstream readFile(toCString(pathToFile));
    sampLenF data;
    int i = 0;
    while(readFile >> data) {
        
        //convert the strings from file into int/double
        int start = atoi(data.start.c_str());
        int end = atoi(data.end.c_str());
        double mean = stod(data.mean);
        double sd = stod(data.sd);
         
        appendValue (grp, start);
        appendValue (grp, end);
        appendValue (grp, mean);
        appendValue (grp, sd);
        appendValue (sampLenPer, grp);
        clear (grp);
    }
    
    return sampLenPer;
}


void readSeqFromFastaFile ( String<char> const & pathToFile, StringSet<String<AminoAcid> > & seqSet, StringSet<String<char> > & seqIds)
{

    //String<char> seqFileName = getAbsolutePath("demos/tutorial/sequence_io/example.fa");
    String<char> seqFileName = getAbsolutePath ( toCString ( pathToFile ) );
    

    SeqFileIn seqFileIn;
    if ( !open ( seqFileIn, toCString ( seqFileName ) ) ) {
        std::cerr << "ERROR: Could not open the file.\n";
        return;
    }

    try {
        readRecords ( seqIds, seqSet, seqFileIn );

    } catch ( Exception const & e ) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return;
    }
    
    /*
    for ( unsigned i = 0; i < length ( seqIds ); ++i ) {
        std::cout << seqIds[i] << '\t' << seqSet[i] << '\n';
    }*/

    //appendValue(seqSet, pathToFile);


    return ;
}

void writeRegsIntoFile (String<String<String<int> > > const & segments, StringSet<String<char> > seqIds, String<char> const & pathToFile) {
    String<char> seqFileName = getAbsolutePath ( toCString ( pathToFile ) );
    
    /*SeqFileOut seqFileOut;
    if ( !open ( seqFileOut, toCString ( seqFileName ) ) ) {
        std::cerr << "ERROR: Could not open the file.\n";
        return ;
    }*/
    
    ofstream seqFileOut;
    seqFileOut.open (toCString(pathToFile));
    
    try {
        for (unsigned i = 0; i < length(segments); ++i) {
            seqFileOut << ">" <<toCString(seqIds[i]) << std::endl;
            for (unsigned x = 0; x < length(segments[i]); ++x) {
                seqFileOut << segments[i][x][0] << " - " << segments[i][x][1] << std::endl;
            }
            
        }
    } catch (Exception const & e ) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return;
    }
    seqFileOut.close();
    return;
}
 
  
/* 2-gram complexity measure
 * the entire sequence is considered as a permutation of 2-grams
 */
template <size_t rows, size_t cols>
void computeTwoLetterScoringMatrix (double (&array)[rows][cols], StringSet<String<char> > const & twoGramAlphabet)
{
    
    double twoLetScore;
    for (unsigned i = 0; i < 676; ++i) {
        for (unsigned j = 0; j < 676; ++j) {
            twoLetScore = 0;
            String<char> lets1 = twoGramAlphabet[i];
            String<char> lets2 = twoGramAlphabet[j];
            
            twoLetScore = twoLetScore + score(blosum62Matrix, lets1[0], lets2[0]);
            twoLetScore = twoLetScore + score(blosum62Matrix, lets1[1], lets2[1]);
            
            array[i][j] = twoLetScore / 2;
        }
    }
    
    return;
}


int strToInt (string str) {
    int first = str[0];
    int second = str[1];
    int trans = (first-65)*26+(second-65);
    
    return trans;
}


template <size_t rows, size_t cols>
double computeTwoLetterEntropy(String<AminoAcid> const & subString, double (&twoLetterScoringMatrix)[rows][cols]) {
    
    double entropy = 0.;
    //TODO useless: ?!StringSet<String<char> > twoGramAlphabet = initTwoGramAlphabet();
    String<int> twoLetterFrequencies;
    resize (twoLetterFrequencies, 676);
    
    
    for (unsigned i = 0; i < 676; ++i) {
        twoLetterFrequencies[i] = 0;
    }
    
    for (unsigned i = 0; i < length(subString)-1; ++i) {
        String<char> twoLet = subString[i];
        appendValue(twoLet, subString[i+1]);
        ++twoLetterFrequencies[strToInt(toCString(twoLet))];
    }
    
    // total number of different 2-letter combinations in seq
    int k = 0; 
    double singleProb = 0;
    double sum = 0;
    String<double> ps;
    resize (ps, 676);
    //initialize ps with 0's
    for (unsigned i = 0; i < 676; ++i) {
        ps[i] = 0.;
    }
    
    for (unsigned i = 0; i < 676; ++i) {
        singleProb = 0;
        if (twoLetterFrequencies[i] != 0) {
            for (unsigned j = 0; j < 676; ++j) {
                if (twoLetterFrequencies[j] != 0) {
                    
                    singleProb = singleProb + (twoLetterFrequencies[j] * twoLetterScoringMatrix[i][j]);
                }
            }
            ps[k] = singleProb;
            sum = sum + singleProb;
            k += 1;
        }
    }
    
    
    for (unsigned i = 0; i < k; ++i) {
        ps[i] = ps[i] / sum;
        entropy = entropy + ( 0. - ps[i] * log( ps[i] ));   
    }
    
    
    //std::cout << "ENTROPY: ";
    //printSeqAnString(subString);
    //std::cout << ": " << entropy << std::endl;
    
    return entropy;
}

/* Extend current blocks to left or right in single steps 
 * if the extended region has a lower entropy, we replace the old region
 */  
String<int> extend (int startPos, int endPos, int limit, string const direction, String<AminoAcid> const & seq) {
    double comp1 = 0;
    double comp2 = 0;
    int pointer = 0;
    int startDecPos = startPos;
    int endDecPos = endPos;
    
    String<AminoAcid> extReg = infix(seq, startPos-1, endPos);
    String<int> decRegs;
    
    if (direction == "left") {
        
        bool dec = false;
        pointer = startPos -2;
        while ((pointer > limit) && (pointer > (startPos - 17))) {
            
            comp1 = computeTwoLetterEntropy(extReg, twoLetterScoringMatrix);
            extReg = infix(seq, pointer, endPos);
            comp2 = computeTwoLetterEntropy(extReg, twoLetterScoringMatrix);
            
            if (comp1 > comp2) {
                if (!dec) {
                    dec = true;
                }    
            } else if (dec) {
                dec = false;
                startDecPos = pointer + 2;
                if (comp1 < comCut) {
                    appendValue(decRegs, startDecPos);
                    appendValue(decRegs, endDecPos);
                }
            }
            
            pointer = pointer -1;
        }
        
        if ((dec) && (pointer == (startPos - 17))) {
            while ((pointer > limit) && (dec)) {
                
                comp1 = computeTwoLetterEntropy(extReg, twoLetterScoringMatrix);
                extReg = infix(seq, pointer, endPos);
                comp2 = computeTwoLetterEntropy(extReg, twoLetterScoringMatrix);
                                            
                if (comp1 < comp2) {
                    dec = false;
                    startDecPos = pointer + 2;
                    if (comp1 < comCut) {
                        appendValue(decRegs, startDecPos);
                        appendValue(decRegs, endDecPos);
                    }
                }
                --pointer;
            }
        }
    
    // the left extension touches the end of the last block of the
    // current lcr blocks
        if ((pointer == limit) && (dec)) {
            startDecPos = pointer + 2;
            if (comp2 < comCut) {
                appendValue(decRegs, startDecPos);
                appendValue(decRegs, endDecPos);
            }
        }
        if (length(decRegs) == 0) {
            appendValue(decRegs, startPos);
            appendValue(decRegs, endPos);
            
        } else {
            //std::cout << "left: NOT empty" << std::endl;
            
        }
    
    } else {
                
        bool dec = false;
        pointer = endPos+1;
        
        while ((pointer < limit) && (pointer < (endPos + 15))) {
            
            comp1 = computeTwoLetterEntropy(extReg, twoLetterScoringMatrix);
            extReg = infix(seq, startPos - 1 , pointer);
            comp2 = computeTwoLetterEntropy(extReg, twoLetterScoringMatrix);
            
            if (comp1 > comp2) {
                if (!dec) {
                    dec = true;
                }    
                
            } else if (dec) {
                dec = false;
                endDecPos = pointer -1;
                if (comp1 < comCut) {
                    appendValue(decRegs, startDecPos);
                    appendValue(decRegs, endDecPos);
                }
            }
            
            ++pointer;
        }
       
        if ((dec) && (pointer == (endPos + 15))) {
            while ((pointer < limit) && (dec)) {
                
                comp1 = computeTwoLetterEntropy(extReg, twoLetterScoringMatrix);
                extReg = infix(seq, startPos - 1, pointer);
                comp2 = computeTwoLetterEntropy(extReg, twoLetterScoringMatrix);
                                            
                if (comp1 < comp2) {
                    dec = false;
                    endDecPos = pointer - 1;
                    if (comp1 < comCut) {
                        appendValue(decRegs, startDecPos);
                        appendValue(decRegs, endDecPos);
                    }
                }
                
                ++pointer;
            }
            
            if ((pointer == limit) && (dec)) {
                endDecPos = limit - 1;
                appendValue(decRegs, startDecPos);
                appendValue(decRegs, endDecPos);
            }
        }
        
        if (length(decRegs) == 0) {
            appendValue(decRegs, startPos);
            appendValue(decRegs, endPos);
        } else {
            //std::cout << "right: NOT empty" << std::endl;
        }        
    }  
    
    return decRegs;
}

/* checks if an element exists
 * 
 */
bool stringContainsItem (String<int> const & str, int id) {
    for (unsigned i = 0; i < length(str); ++i) {
        if (str[i] == id) {
            return true;
        }
    }
    return false;
}

/* 
 * 
 */
void convertIdsToPosition(String<int> & lcrBlocks, String<int> const & uniquePathVertIds, String<String<int> > const & vertexArray) {
    for (unsigned i = 0; i < length(uniquePathVertIds); ++i) {
        if (!stringContainsItem(lcrBlocks, vertexArray[uniquePathVertIds[i]][1]) && vertexArray[uniquePathVertIds[i]][1] != 0) {
            appendValue(lcrBlocks, vertexArray[uniquePathVertIds[i]][1]);
            
        }
        if (!stringContainsItem(lcrBlocks, vertexArray[uniquePathVertIds[i]][2]) && vertexArray[uniquePathVertIds[i]][2] != 0) {
            appendValue(lcrBlocks, vertexArray[uniquePathVertIds[i]][2]);
            
        }
        
    }
    
    int len = length(lcrBlocks);
    quickSort(lcrBlocks, 0, len-1);
    
    return;
}

/* returns regions blockwise by position in the sequence e.g. 4-16
 * 
 */
String<String<int> > genBlocksFromList(String<int> const & blockList) {
    String<int> tmpBlockList = blockList;
    String<String<int> > lcrBlocks;
    int pos = 0;
    int counter = 0;
    int startPosLcrBlock = tmpBlockList[pos];    
    
    std::cout << "liste" << std::endl;
    printSeqAnString(tmpBlockList);
    
    for (unsigned i = 0; i < length(tmpBlockList) ; ++i) {
        if (tmpBlockList[i]+1 !=tmpBlockList[i+1] ) {
            int endPosLcrBlock = tmpBlockList[i];
            appendValue(lcrBlocks, startPosLcrBlock);
            appendValue(lcrBlocks[counter], endPosLcrBlock);
            //std::cout << startPosLcrBlock << "-" << endPosLcrBlock << std::endl;
            startPosLcrBlock = tmpBlockList[pos+1];
            ++counter;
        }
        ++pos;
    }
    
    std::cout << "lcrBlocks: " << std::endl;
    for (unsigned i = 0; i < length(lcrBlocks); ++i) {
        printSeqAnString(lcrBlocks[i]);
    }
    
    return lcrBlocks;
}

/* 
 * 
 */
String<String<int> > generateLCRBlocks (String<int> const & longPaths, String<String<int> > const & vertexArray) {
    String<int> lcrBlockList;
    String<int> uniquePathVertIds;
    
    int pathLength = length(longPaths);
    for (unsigned i = 0; i < pathLength; ++i) {
        if (!stringContainsItem(uniquePathVertIds, longPaths[i])) {
            appendValue (uniquePathVertIds, longPaths[i]);
        }
    }
    
    convertIdsToPosition(lcrBlockList, uniquePathVertIds, vertexArray);
    String<String<int> > lcrBlocks = genBlocksFromList (lcrBlockList);
    
    return lcrBlocks;
}

/* working on LCR blocks (extending while considering limits to the left and right)
 * limits are the end of the sequence or found blocks positions
 */
String<String<int> > processBlocks (String<String<int> > const & blocks, String<AminoAcid> const & seq) {
    String<int> tmpBlocks;
    resize (tmpBlocks, 2);
    String<int> saveExtendedBlock;
    String<String<int> > extendedBlocks;
    int limit = 0;
    int lcrBlockStart = 0;
    int lcrBlockEnd = 0;
    
    
    for (unsigned i = 0; i < length(blocks); ++i) {
        if (i == 0) {
            limit = -1;
        } else {
            limit = lcrBlockEnd -1;
        }
        
        tmpBlocks = extend (blocks[i][0], blocks[i][1], limit, "left", seq);
        appendValue (saveExtendedBlock, tmpBlocks[0]);
        
        limit = length(seq);
        
        tmpBlocks = extend (blocks[i][0], blocks[i][1], limit, "right", seq);
        lcrBlockEnd = tmpBlocks[1];
        
        appendValue (saveExtendedBlock, tmpBlocks[1]);
        appendValue(extendedBlocks, saveExtendedBlock);
        clear(tmpBlocks);
        clear(saveExtendedBlock);
    }
    
    return extendedBlocks;
}

/* merging overlapping regions e.g. 3-14 and 12-22 -> 3-22
 * 
 */
String<String<int> > mergeRegions (String<String<int> > const & extendedBlocks) {
    String<String<int> > lcrs;
    String<int> initLcrs;
    appendValue(initLcrs, extendedBlocks[0][0]);
    appendValue(initLcrs, extendedBlocks[0][1]);
    appendValue(lcrs, initLcrs);
    String<int> tmpBlock;
    int counter = 0;
    
    for (unsigned i = 0; i < length(extendedBlocks)-1; ++i) {
        int curBlockStart = lcrs[counter][0];
        int curBlockEnd = lcrs[counter][1];
        
        if (i == 0) {
            int curBlockStart = extendedBlocks[0][0];
            int curBlockEnd = extendedBlocks[0][1];
        }
        
        int nextBlockStart = extendedBlocks[i+1][0];
        int nextBlockEnd = extendedBlocks[i+1][1];
        
        if ((curBlockStart <= nextBlockStart) && (nextBlockStart <= curBlockEnd)) {
            int lcrStart = min(curBlockStart, nextBlockStart);
            int lcrEnd = max(curBlockEnd, nextBlockEnd);
            
            lcrs[counter][0] = lcrStart;
            lcrs[counter][1] = lcrEnd;
            
        } else {
            int lcrStart = nextBlockStart;
            int lcrEnd = nextBlockEnd;
            appendValue (tmpBlock, lcrStart);
            appendValue (tmpBlock, lcrEnd);
            appendValue (lcrs, tmpBlock);
            clear (tmpBlock);
            ++counter;
        }
    }
    
    return lcrs;
}

/* checkLeftRegs() is checking regions on the left side of the calculated alignment and adds regions, whose complexity
 * is lower than a calculated cutoff value(cCut)
 * 
 */
StringSet<String<int> > checkLeftRegs (int const & aliStart, int const & aliEnd, int const & start, int const & end, String<AminoAcid> const & seq, double const & cCut) {
    double com = 0;
    StringSet<String<int> > left;
    String<int> leftPrep;
        
    if (aliStart > 7) {     
        String<AminoAcid> str = infix(seq, start-1, start+aliStart-2);
        
        com = computeTwoLetterEntropy(str, twoLetterScoringMatrix);
        com = com/(length(str)-1);
        
        if (com <= cCut) {
            appendValue (leftPrep, start);
            appendValue (leftPrep, (start+aliStart-2));
            appendValue (left, leftPrep);
            clear (leftPrep);
        }    
    }
    
    int ergebnis = end-start+1-aliEnd;
    
    
    if (ergebnis > 7) {
        
        String<AminoAcid> str = infix(seq, start+aliEnd-1, end);
        com = computeTwoLetterEntropy(str, twoLetterScoringMatrix);
        com = com/(length(str)-1);
        if (com <= cCut) {
            appendValue (leftPrep, start+aliEnd);
            appendValue (leftPrep, end);
            appendValue (left, leftPrep);
            
            clear (leftPrep);
        }
    }
    
    return left;
}

/* addToResult() adds calculated left-regions to the StringSet result
 * 
 */
StringSet<String<int> > addToResult (StringSet<String<int> > & result, StringSet<String<int> > const & left) {
    int j = 0;
    
    for (unsigned i = 0; i <length(left); ++i) {
        String<int> str1 = left[i];
        int endLeft = str1[1];
        bool found = false;
        while (!found) {
            if (j < length(result)) {
                String<int> str2 = result[j];
                int startResult = str2[0];
                if (endLeft < startResult) {
                    found = true;
                    insertValue (result, j, str1);
                    j = j + 2;
                } else {
                    ++j;
                }
            } else {
                appendValue (result, str1);
                
                found = true;
            }
        }
    }
    
    return result;
}

/* findAlignment() finds a local alignment (scoringmatrix: Blosum62) with the Smith-Waterman algorithm
 * and returns the position of the query found in the sequence
 * format of aliPos: start1 end1 start2 end2
 */
StringSet<String<int> > findAlignment (String<AminoAcid> const & seq1, String<AminoAcid> const & seq2) {
    typedef String<AminoAcid> TSequence;
    typedef Align<TSequence> TAlign;
    
    String<int> pos;
    StringSet<String<int> > aliPos;
    
    TSequence seqH = seq1;
    TSequence seqV = seq2;
    
    try {
        TAlign align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);
        
        //matrixScore, gapcosts: open: 10, extension:0.5
        int result = localAlignment(align, Score<int, ScoreMatrix<AminoAcid, ScoreSpecBlosum62> >(-0.5,-10));
        int len=length(row(align,0));
        
        
        /*
        std::cout << "            SequenPos: 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7" << std::endl;
        std::cout << "            Sequence1: ";
        printSeqAnString(seq1) ;
        std::cout << "            Sequence2: ";
        printSeqAnString(seq2);
        std::cout << "            Score: " << result << "\n";
        std::cout << "            The resulting alignment is\n" << align << "\n";
        */
        
        auto src1Pos = toSourcePosition(row(align,0),0);
        auto src1End = toSourcePosition(row(align,0),len-1);
        auto src2Pos = toSourcePosition(row(align,1),0);
        auto src2End = toSourcePosition(row(align,1),len-1);
        
        int similarity = 0;
        
        
        for (unsigned i = 0; i < len; ++i) {
            AminoAcid firstChar = row(align,0)[i];
            AminoAcid secChar = row(align,1)[i];
            
            //print score of aligned characters
            //std::cout << "            " <<(char)firstChar << "-" << (char)secChar << " has a score of: " << score(blosum62Matrix, firstChar, secChar) << std::endl;
            
            if (firstChar == secChar) {
                ++similarity;
            } else if (score(blosum62Matrix, firstChar, secChar) > 0) {
                ++similarity;
            }
        }
        
        
        if (similarity > 4) {
            appendValue (pos, src1Pos);
            appendValue (pos, src1End);
            appendValue (aliPos, pos);
            clear (pos);
            appendValue (pos, src2Pos);
            appendValue (pos, src2End);
            appendValue (aliPos, pos);
        }
    } catch (Exception const & e) {
        std::cout << "ERROR: " << e.what() << std::endl;
    }
    
    return aliPos;
}

/*
 * mark == -1 means forward, mark == -2 means backwards 
 */
StringSet<String<int> > checkAdjBlock (int const & start1, int const & end1, String<int> adjBlock, String<AminoAcid> const & seq, double const & cCut, int const & mark) {    
    double com = 0;
    StringSet<String<int> > result;
    
    //TODO better ? infix(seq, start1-2, ...);
    String<AminoAcid> seq1 = infix(seq, start1-1, end1);
    int start2 = adjBlock[0];
    int end2 = adjBlock[1];
    //TODO better?  infix(seq, start2-2, ...);
    String<AminoAcid> seq2 = infix(seq, start2-1, end2);
    
    StringSet<String<int> > aliPos;
    if (mark == -1) {
        //TODO where'S the difference between seq2,seq1 and seq1,seq2 ??!?!
        aliPos = findAlignment (seq2, seq1);
    } else {
        aliPos = findAlignment (seq1, seq2);
    }
    
    /*
    for (unsigned i = 0; i < length(aliPos); ++i) {
        std::cout << "        ";
        printSeqAnString(aliPos[i]);
    }*/
    
    if (length(aliPos) != 0) {
        int aliStart2 = aliPos[0][0];
        int aliEnd2 = aliPos[0][1];
            
        int aliStart1 = aliPos[1][0];
        int aliEnd1 = aliPos[1][1];
        //TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO         
        
        int append1 = 0;
        int append2 = 0;
        
        if (mark == -1) {
            
            append1 = start1 + aliStart1;
            append2 = start1 + aliEnd1;
            
            String<int> resultPrep;
            appendValue (resultPrep, append1);
            appendValue (resultPrep, append2);
            
            appendValue (result, resultPrep);
            StringSet<String<int> > left = checkLeftRegs (aliStart1, aliEnd1, start1, end1, seq, cCut);
            /*
            for (unsigned i = 0; i < length(left); ++i) {
                std::cout << "        ";
                printSeqAnString(left[i]);
            }*/
            addToResult (result, left);
            
        } else {
            String<int> resultPrep;
            append1 = start1 + aliStart2;
            append2 = start1 + aliEnd2;
            appendValue (resultPrep, append1);
            appendValue (resultPrep, append2);
            appendValue (result, resultPrep);
            clear (resultPrep);
            
            StringSet<String<int> > left = checkLeftRegs (aliStart2, aliEnd2, start1, end1, seq, cCut);
            /*
            for (unsigned i = 0; i < length(left); ++i) {
                std::cout << "        ";
                printSeqAnString(left[i]);
            }*/
            
            addToResult (result, left);
            
            append1 = start2 + aliStart1 ;
            append2 = start2 + aliEnd1 ;
            appendValue (resultPrep, append1);
            appendValue (resultPrep, append2);
            appendValue (result, resultPrep);
            clear (resultPrep);
            
            left = checkLeftRegs (aliStart1, aliEnd1, start2, end2, seq, cCut);
            /*
            for (unsigned i = 0; i < length(left); ++i) {
                std::cout << "        ";
                printSeqAnString(left[i]);
            }*/
            
            addToResult (result, left);
            
        }
        
    }
    
    return result;
}


/* cheks if there are deletable characters in a region 
 * 
 */
StringSet<String<int> > checkDeletability (String<String<int> > const & lcrs, int const & maxIndex, String<AminoAcid> const & seq, double const & cCut) {    
    StringSet<String<int> > result;
    String<int> block;
    int start1 = 0;
    int end1 = 0;
    
    block = lcrs[maxIndex];
    start1 = block[0];
    end1 = block[1];
    
    if (maxIndex != 0) {
        block = lcrs[maxIndex-1];
        try {
            result = checkAdjBlock (start1, end1, block, seq, cCut, -1);
            
        } catch ( Exception const & e ) {
            std::cout << "ERROR: " << e.what() << std::endl;
            return result;
        }
    }
    
    if (length(result) == 0) {
        if (maxIndex != length(lcrs)-1) {
            block = lcrs[maxIndex+1];
            
            try {
                result = checkAdjBlock (start1, end1, block, seq, cCut, -2);
            } catch ( Exception const & e ) {
               std::cout << "ERROR: " << e.what() << std::endl;
                return result;
            }
            
            if (length(result) != 0) {
                //appending a -2 to result for 'back'
                appendValue (result, -2);
                
            }
        }
    }
    
    return result;
}


/*
 * 
 */
String<String<int> > filter (String<String<int> > const & mergedBlocks, String<AminoAcid> const & seq) {
    
    int seqLen = length(seq);
    int len = length(mergedBlocks);
    String<double> complexityList ;
    String<String<int> > tmpBlocks = mergedBlocks;
    String<String<double > > sampLenPer = readSampleLenPercentage("sampledLenRepPer");
    double comp = 0.;
    double max = -222222222;
    int it = 0;
    
    if (len != 1) {
        for (unsigned i = 0; i < len; ++i) {
            int start = tmpBlocks[i][0];
            int end = tmpBlocks[i][1];
            String<AminoAcid> str = infix(seq, start-1, end);
            
            comp = computeTwoLetterEntropy(str, twoLetterScoringMatrix);
            comp = comp/(length(str)-1);
            appendValue (complexityList, comp);
        }
        
        bool found = false;
        int lineNr = 0;
        double limit=0;
        double preLongest = 0;
        double prePer = 0;
        int maxIndex = 0;
        double cCut = 0;
        double per = 0;
        int shortest = 0;
        int longest = 0;
        
        while (!found && (lineNr <length(sampLenPer))) {
            shortest = sampLenPer[lineNr][0];
            longest = sampLenPer[lineNr][1];
            per = sampLenPer[lineNr][2];
            
            if ((seqLen >= shortest) && (seqLen <= longest)) {
                limit = len * (1-per);
                found = true;
            } else if (seqLen < shortest) {
                int diff1 = shortest - seqLen;
                int diff2 = seqLen - preLongest;
                if ((diff2 > diff1) || (preLongest == 0)) {
                    limit = len * (1-per);
                    found = true;
                } else {
                    limit = len * (1-prePer);
                    found = true;
                }
            }
            
            preLongest = longest;
            prePer = per;
            ++lineNr;
        }
     
        
        if (!found) {
            limit = len * (1 - per);
        }
        lineNr = 0;
        
        while (lineNr < limit) {
            int compLen = length(complexityList);
            
            while (it < compLen) {
                double comp = complexityList[it];
                if (comp > max) {
                    max = comp;
                    maxIndex = it;
                }
                ++it;
            }
            
            cCut = complexityList[maxIndex];
            ++lineNr;
        }
        lineNr = 0;
        it = 0;
        len = length(tmpBlocks);
        StringSet<String<int> > result;
        
        while ((lineNr < limit) && (len != 1) && (it < len)) {
            int start = tmpBlocks[it][0];
            int end = tmpBlocks[it][1];
            String<AminoAcid> str = infix(seq, start-1, end);
            
            comp = computeTwoLetterEntropy(str, twoLetterScoringMatrix);
            comp = comp/(length(str)-1);
            
            if (comp >= cCut ) {
                result = checkDeletability (tmpBlocks, it, seq, cCut);
                
                int rSize = length(result);
                std::cout << "RSIZE: " << rSize << std::endl;
                bool fromBack = false;
                if (rSize != 0) {
                    
                    int lastElem = result[rSize-1][0];
                    // a -2 indicates 'back', -1 indicates 'front'
                    if (lastElem == -2) {
                        --rSize;
                        fromBack = true;
                    }
                    
                    erase(tmpBlocks, it);
                    
                    for (unsigned k = 0; k < rSize; ++k) {
                        String<int> strRange = result[k];
                        insertValue (tmpBlocks, it, strRange);
                        ++it;
                    }
                    
                    if (fromBack) {    
                        erase (tmpBlocks, it);
                    }
                    len = length(tmpBlocks);
                
                } else {
                    
                    erase (tmpBlocks, it);
                    ++lineNr;
                    len = length(tmpBlocks);
                }
            } else {
                ++it;
            }
        }
    }
    
    /*
    std::cout << "MEINE TMPBLOCKS: " << std::endl;
    
    for (unsigned i = 0; i < length(tmpBlocks); ++i) {
        for (unsigned j = 0; j < length(tmpBlocks[i]); ++j) {
            std::cout << tmpBlocks[i][j];
        }
        std::cout << std::endl;
    }*/
    
    
    
    return tmpBlocks;
        
}
    

void segSeq (String<AminoAcid> const & seq, int offset, int const & winLen, float const loCut, float const hiCut, StringSet<String<int> > & segs) {
    
    String<int> seg;
    int downset = (winLen+1)/2 - 1;
    int upset = winLen - downset;
    int first = downset;
    int last = length(seq) - upset;
    int lowLim = first;
    
    
    String<long double> H = initSeqEnt(seq);
    
    if (length(H) == 0) {
        return ;
    }
    
    
    // save calculated shanEntropy in vector H
    for (unsigned i = 0; i <= last; ++i) {
        int winStartPos = i;
        int winEndPos = i+winLen;
        
        if (winEndPos > length(seq)) {
            break;
        }
        
        String<int> freqVector = initFreqVector();
        
        String<AminoAcid> window = infix (seq, winStartPos, winEndPos);
        
        getAminoAcidCountInSequence(freqVector, window);
        
        long double K2KV = shanEntropy (winLen, freqVector);
        
        
        H[i+first] = K2KV;
        
        
    }
   
    //sliding window
    for (int i = first; i <= last; ++i) {
        if (H[i] <= loCut && H[i] != -1) {
            int loi = findLow (hiCut, i, lowLim, H);
            int hii = findHigh (hiCut, i, last, H);
            
            int leftEnd = loi - downset;
            int rightEnd = hii + upset -1;
            
            String<AminoAcid> tempSeq = infix(seq, leftEnd, rightEnd + 1);
            
            trim(tempSeq, leftEnd, rightEnd, winLen);
            
            // check triggered window left side
            if (i+upset-1 < leftEnd) {
                int lend = loi - downset;
                int rend = leftEnd;
                
                String<AminoAcid> leftSeq = infix(seq, lend, rend);
                StringSet<String<int> > leftSegs;
                
                segSeq (leftSeq, offset+lend, winLen, loCut, hiCut, leftSegs);
                
                if (leftSegs != 0) {
                    appendValue (segs, leftSegs[0]);
                }
            }
            
            String<int> region;
            appendValue (region, leftEnd+offset+1);
            appendValue (region, rightEnd+offset+1);
            appendValue (segs, region);
            
            i = std::min(hii, rightEnd+downset);
            lowLim = i + 1;
            
        }
    }
    
    return ;
}



// ---------------------------          SEG
String<String<int> > seg ( String<AminoAcid> const & seq, unsigned const & winLen, float const & K2A, float const & K2B )
{
    int offset = 0;
    
    StringSet<String<int> > segs;
    
    segSeq (seq, offset, winLen, K2A, K2B, segs);
    
    if (length(segs)==0) {
        std::cout << "No low-complexity-Regions found" << std::endl;
        return segs;
    } 
    
    String<String<int> > mergedRegions = mergeRegions(segs);
    std::cout << "low complexity regions found at: " << std::endl;
    for (unsigned i = 0; i < length(mergedRegions); ++i) {
        std::cout << mergedRegions[i][0] << " - " << mergedRegions[i][1] << '\t';
        /*
        for (int x = mergedRegions[i][0]-1; x < mergedRegions[i][1]; ++x) {
         
            std::cout << seq[x];
        }*/
        std::cout << std::endl;
    }
    
    
    return mergedRegions;
}

// ---------------------------          GBA
int gba (String<AminoAcid> const & seq, int const winLen, int const maxIndelTandem, int const maxIndelCryptic) {
    
    String<int> vertex;
    String<String<int> > vertexArray;

    if ( winLen > length ( seq ) | maxIndelTandem > winLen | maxIndelCryptic > winLen ) {
        std::cerr << "Input-parameters are not allowed." << std::endl;
        return 0;
    }
    
    TGraph gbaGraph = fastCreateVertices ( seq, winLen, maxIndelTandem, maxIndelCryptic, vertex, vertexArray );
    gbaGraph = createEdges ( gbaGraph, vertexArray, seq );
    //printVertexArray ( vertexArray );
    
    String<int> longPathVertices = getLongestPaths ( gbaGraph );

    int rows = 27;
    int columns = 27;

    String<int> wholeSeqFrequencies = initFreqVector();

    
    getAminoAcidCountInSequence(wholeSeqFrequencies, seq);

    String<String<int> > lcrBlocks = generateLCRBlocks (longPathVertices, vertexArray);
    
    
    String<String<int> > extendedBlocks = processBlocks(lcrBlocks, seq);

    
    String<String<int> > mergedBlocks = mergeRegions(extendedBlocks);
    
    String<String<int> > filteredRegions = filter(mergedBlocks, seq);
    filteredRegions = filter(filteredRegions, seq);
    filteredRegions = mergeRegions(filteredRegions);
    
    
    std::cout << "low complexity regions found at: " << std::endl;
    
    for (unsigned i = 0; i < length(filteredRegions); ++i) {
        std::cout << filteredRegions[i][0] << " - " << filteredRegions[i][1] << std::endl;;
    }
    
    return 1;
    
    
}
 
struct Options
{
    seqan::CharString inputFile;
    seqan::CharString outputFile;
    seqan::CharString method;
    int window;
    float locut;
    float hicut;
    int maxIndelTandem;
    int maxIndelCryptic;
};
 

 

seqan::ArgumentParser::ParseResult parseCommandLine(Options & options, int argc, char const ** argv) {
    seqan::ArgumentParser parser("find_lcr");
    
    setShortDescription(parser, "Identifying low-complexity Regions");
    addDescription(parser,
                   "The program will identify low-complexity Regions in proteinsequences using GBA or SEG");
    addArgument(parser, seqan::ArgParseOption("m", "method", "Method to identify LCRs: GBA, SEG", seqan::ArgParseOption::STRING, "METHOD"));
    addArgument(parser, seqan::ArgParseOption("s", "seq", "FASTA file with one or more sequences.", seqan::ArgParseOption::INPUT_FILE, "INPUT"));
    addArgument(parser, seqan::ArgParseOption("s", "out", ".txt file with low-complexity regions.", seqan::ArgParseOption::OUTPUT_FILE, "OUTPUT"));
    addOption(parser, seqan::ArgParseOption("w", "window", "Length of sliding window.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "window", 12);
    addOption(parser, seqan::ArgParseOption("lc", "locut", "Locut parameter for SEG.", seqan::ArgParseArgument::DOUBLE, "FLOAT"));
    setDefaultValue(parser, "locut", 2.2);
    addOption(parser, seqan::ArgParseOption("hc", "hicut", "Hicut parameter for SEG.", seqan::ArgParseArgument::DOUBLE, "FLOAT"));
    setDefaultValue(parser, "hicut", 2.5);
    
    addOption(parser, seqan::ArgParseOption("tr", "indelTandemRepeats", "Parameter to specify the maximum number of insertions and deletions between similar repeats for GBA.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "indelTandemRepeats", 3);
    addOption(parser, seqan::ArgParseOption("cr", "indelCrypticRepeats", "Parameter to specify the maximum distance between letters in cryptic repeats for GBA.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "indelCrypticRepeats", 5);
    
    
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    
    if (res != seqan::ArgumentParser::PARSE_OK) {
        return res;
    }
    
    getArgumentValue(options.method, parser, 0);
    getArgumentValue(options.inputFile, parser, 1);
    getArgumentValue(options.outputFile, parser, 2);
    getOptionValue(options.window, parser, "window");
    getOptionValue(options.locut, parser, "locut");
    getOptionValue(options.hicut, parser, "hicut");
    getOptionValue(options.maxIndelTandem, parser, "indelTandemRepeats");
    getOptionValue(options.maxIndelCryptic, parser, "indelCrypticRepeats");
    
    return seqan::ArgumentParser::PARSE_OK;
}  
  

int main ( int argc, char const ** argv )
{
    
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    
    
    String<String<int> > segs;
    String<String<String<int> > > segments;
    
    StringSet<String<char> > seqIds;
    StringSet<String<AminoAcid> > seqSet;
    String<AminoAcid> seqPaper = "GAYTSVAYTVPQAWTVW"; //FFFFFFFFPNGKLJHAMFFFFFQE";
    StringSet<String<char> > twoGramAlphabet = initTwoGramAlphabet();
    
    
    String<char> alg = options.method;
    String<char> fileIn = options.inputFile;
    String<char> fileOut = options.outputFile;
    int windowLength = options.window;
    float loCut = options.locut;
    float hiCut = options.hicut;
    int maxIndelTandem = options.maxIndelTandem;
    int maxIndelCryptic = options.maxIndelCryptic;
    
    
    computeTwoLetterScoringMatrix(twoLetterScoringMatrix, twoGramAlphabet);
    normalizeScoringMatrix(twoLetterScoringMatrix);
    
    if (alg == "gba" or alg == "GBA") {
        readSeqFromFastaFile ( fileIn, seqSet, seqIds );
        std::cout << "------------------------------- executing GBA... -----------------------------" << std::endl;
        SEQAN_OMP_PRAGMA (parallel for)
        for (unsigned i = 0; i < length(seqSet); ++i) { 
            String<AminoAcid> seq = seqSet[i];
            
            gba(seq, windowLength, maxIndelTandem, maxIndelCryptic);
            appendValue (segments, segs);
            /*
            std::cout << "1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0" << std::endl;
            printSeqAnString(seqPaper);
            */
        }
        
    } else if (alg == "seg" or alg == "SEG") {
        readSeqFromFastaFile ( fileIn, seqSet, seqIds );
        std::cout << "------------------------------- executing SEG... -----------------------------" << std::endl;
        SEQAN_OMP_PRAGMA (parallel for)
        for (unsigned i = 0; i < length(seqSet); ++i) {
            String<AminoAcid> seq = seqSet[i];
            segs = seg(seq, windowLength, loCut, hiCut);
            appendValue (segments, segs);
        }
        
    } else {
        std::cout << "Please specify whether you want to run 'gba' or 'seg'" << std::endl;
        return 1;
    }
    
    writeRegsIntoFile(segments, seqIds, fileOut);

    
    
    
    
    //high_resolution_clock::time_point t1 = high_resolution_clock::now();
    //do stuff
    //high_resolution_clock::time_point t2 = high_resolution_clock::now();

    //std::cout << "gemessene Dauer mit functuon: " << funcTime(createVertices, seq, winLen, maxIndelTandem, maxIndelCryptic) << std::endl;

    
   //for (unsigned i = 0; i < length(extendedBlocks); ++i) {       
       //maskSequence(seq, extendedBlocks[i][0], extendedBlocks[i][1]);
   //}
   

    String<AminoAcid> bla = "GDFBL";
    String<AminoAcid> blaInfix = infix(bla, 2, 4);
    std::cout << blaInfix << std::endl;
    
    
    return 0;
}
