#include <iostream>
#include <fstream>
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

//profiler:
//valgrind --tool=callgrind ./my_proj ...

using namespace std::chrono;
using namespace std;
using namespace seqan;

typedef unsigned int TCargo;
typedef Graph<Directed<TCargo> > TGraph;
typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
typedef Size<TGraph>::Type TSize;
typedef double TScoreValue;

//complexity cutoff value for extending longest path intervals
//the cutoff value was calculated by Li & Kahveci by randomly sampling sequences from Swissprot that contain repeat regions
// (Li & Kahveci, Bioinformatics 2006)
double comCut = 3.7491631225622157;

//alphabet size, 20 AminoAcids, String<AminoAcid> contains 27 letters (wildcards..)
int const alphSize = 20;

//reading the R-NR-Matrices from a file
Score<TScoreValue, ScoreMatrix<AminoAcid, Default> > readRNRMatrix ( String<char> RNRMatrixPath )
{
    Score<TScoreValue, ScoreMatrix<AminoAcid, Default> > RNRMatrix;
    fstream seqFileIn;
     //SeqFileIn seqFileIn;
    seqFileIn.open(toCString ( getAbsolutePath(toCString(RNRMatrixPath))));
    if (seqFileIn.is_open()) {
        seqFileIn.close();
    } else {
        std::cout << "ERROR: Could not open the RNR Matrices from path: " << toCString(getAbsolutePath(toCString(RNRMatrixPath))) << std::endl;
        std::cout << "The vertexconstruction will not work properly and will return a different graph!" << std::endl;
        
    }
    
    /* if ( !open ( seqFileIn, toCString ( getAbsolutePath(toCString(RNRMatrixPath))))) {
         std::cerr << "ERROR: Could not open the RNR Matrices" << "\n";
     }*/

    try {
        loadScoreMatrix ( RNRMatrix, toCString ( getAbsolutePath ( toCString ( RNRMatrixPath ) ) ) );

    } catch ( Exception const & e ) {
        std::cout << "ERROR: " << e.what() << std::endl;
    }

    return RNRMatrix;
}

/* initialize the alphabet vector for calculating the frequency of characters in the sliding window
 * 
 */
String<int> initFreqVector () {
    String<int> alphabetFreq;
    int N=27;
    resize (alphabetFreq, N);
    for (unsigned i = 0; i < N; i++) {
        alphabetFreq[i] = 0;
    }

    return alphabetFreq;
}

/* This is for gba, unneccessary size of the alphabet, do it like initFreqVector (above) but check fastCreateVertices() !!! 
 * 
 */
String<int> initFreqVector2 () {
    String<int> alphabetFreq;
    //sized 123 because the last letter 'z' in ASCII is at position 122
    
    int N=123;
    resize (alphabetFreq, N);
    for (unsigned i = 0; i < N; i++) {
        alphabetFreq[i] = 0;
    }

    return alphabetFreq;
}


//GLOBAL MATRICES
//TODO ADJUST PATHS
String<char> nonRepeatMatrixPath = "../LCR-finder/matrices/combinedMatricesRowByRow095NonRepeat";
String<char> repeatMatrixPath = "../LCR-finder/matrices/combinedMatricesRowByRow095Repeat";

Score<TScoreValue, ScoreMatrix<AminoAcid, Default> > nonRepeatMatrix = readRNRMatrix ( nonRepeatMatrixPath );
Score<TScoreValue, ScoreMatrix<AminoAcid, Default> > repeatMatrix = readRNRMatrix ( repeatMatrixPath );
const Blosum62 blosum62Matrix;
double twoLetterScoringMatrix[676][676] = { {} };
double nonNormalizedTwoLetScoMatrix[676][676] = { {} };
String<int> iniFreqVector = initFreqVector();
String<int> iniFreqVector2 = initFreqVector2();


/* Initiating 2-gram String for any combination of letters in the alphabet
 * 
 */
StringSet<String<char> > initTwoGramAlphabet() {
    String<char> alphabet;
    StringSet<String<char> > alphabetTwoGram;
    resize (alphabet, 27);
    
    for (unsigned i = 0; i < 26; ++i) {
        alphabet[i] = i + 65;
    }
    
    for (unsigned i = 0; i < 26; ++i) {
        //are similar letters also a k-gram?
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
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            if (array[i][j] != 0) {                
                array[i][j] = pow (2, array[i][j]);
            }
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


/* create Edges in our Graph with specific conditions
 * cutoff = 1, tTwo = 3, tThree = 5
 * TODO
 */
TGraph createEdges ( TGraph gbaGraph, String<String<int> > & vertexArray, String<AminoAcid> const & seq ) {
    int cutoff = 1;
    //TODO parameter as input
    int tTwo = 3;
    int tThree = 5;
    int tFour = 15;
    int lenArr = length(vertexArray);

    if ( length ( vertexArray ) <= 1 ) {
        return gbaGraph;
    }

    for ( int startVertex = 0; startVertex < length ( vertexArray ); ++startVertex ) {
        //std::cout << "STARTVERTEX: " << vertexArray[startVertex][1] << "," << vertexArray[startVertex][2] << std::endl;
        for ( int destinationVertex = startVertex + 1; destinationVertex < length ( vertexArray ); ++destinationVertex ) {
            //std::cout << "DESTINATIONVERTEX: " << vertexArray[destinationVertex][1] << "," << vertexArray[destinationVertex][2] << std::endl;
            std::tuple<char, char> first ( seq[vertexArray[startVertex][1]-1], seq[vertexArray[destinationVertex][1]-1] );
            std::tuple<char, char> second ( seq[vertexArray[startVertex][2]-1], seq[vertexArray[destinationVertex][2]-1] );
            
            if ( score ( blosum62Matrix, std::get<0> ( first ), std::get<0> ( second ) ) > cutoff && score ( blosum62Matrix, std::get<1> ( first ), std::get<1> ( second ) ) > cutoff ) {

                // 4 conditions to be checked before creating an edge
                //condition #1:
                if ( abs ( ( vertexArray[startVertex][2] - vertexArray[startVertex][1] ) - ( vertexArray[destinationVertex][2] - vertexArray[destinationVertex][1] ) ) <= tTwo ) {
                    //condition #2:
                    if ((vertexArray[destinationVertex][1] - vertexArray[startVertex][1]) <= tFour) {                    
                    //condition #3:
                        if ( (vertexArray[startVertex][1] <= vertexArray[destinationVertex][1]) && (vertexArray[destinationVertex][1] <= vertexArray[startVertex][2]) && (vertexArray[startVertex][2] <= vertexArray[destinationVertex][2]) ) {
                            //condition #4:
                            if ((vertexArray[startVertex][1] == vertexArray[destinationVertex][1]) || (vertexArray[startVertex][2] == vertexArray[destinationVertex][1]) || (vertexArray[startVertex][2] == vertexArray[destinationVertex][2])) {
                                if ( ((vertexArray[startVertex][2]-vertexArray[startVertex][1]) <= tThree) && (vertexArray[destinationVertex][2]-vertexArray[destinationVertex][1]) <= tThree ) {
                                    addEdge ( gbaGraph, startVertex, destinationVertex );
                                }
                            } else {
                                addEdge ( gbaGraph, startVertex, destinationVertex );
                            }
                        }
                    }
                }
            }
        }
    }
    //std::cout << "number of edges: " << numEdges(gbaGraph) << std::endl;;
    //std::cout << "vertexArray: " << std::endl;
    //for (int i = 0; i < length(vertexArray); ++i) {
    //    printSeqAnString(vertexArray[i]);
    //}
    
    //drawGraph ( gbaGraph );
    return gbaGraph;
}


void createVertexArray ( String<int> & vertex, String<String<int> > & vertexArray, int const counter, int & firstPositionOfTuple, int & secondPositionOfTuple) {
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
bool checkSameWindowByChance ( Score<TScoreValue, ScoreMatrix<TSequenceValue, TSpec> > const & repeatMatrix, Score<TScoreValue, ScoreMatrix<TSequenceValue, TSpec> > const & nonRepeatMatrix, char const & firstLetter, char const & secondLetter, double const & alphabetFreq, double const & winLen ) {
    bool construct = false;
    double scoreRepeatMatrix = score ( repeatMatrix,firstLetter, secondLetter);
    
    double scoreNonRepeatMatrix = score ( nonRepeatMatrix, firstLetter, secondLetter);
    double frequency = alphabetFreq/winLen;

    double calcDiffNR = abs(scoreNonRepeatMatrix - (frequency));
    double calcDiffR = abs(scoreRepeatMatrix - (frequency));
    
    //std::cout << "Score in repeatMatrix: " << scoreRepeatMatrix << " Score nonRepeatMatrix: " << scoreNonRepeatMatrix << " freq: " << frequency << std::endl;
    
    if (calcDiffNR >= calcDiffR) {
        construct = true;
    }
    
    return construct;
}

/* calculate the frequency of each character in the given sequence
 * 
 */
void getAminoAcidCountInSequence (String<int> & wholeSeqFrequencies, String<AminoAcid> const & seq) {
    for (unsigned i = 0; i < length(seq); ++i) {
        char letter = seq[i];
        int idx=(int)letter;
        ++wholeSeqFrequencies[idx-65];
    }
    
    return;
}


//String<AminoAcid> const & seq, String<int> & vertex, StringSet<String<int> > & vertexArray
/* create a vertex for each 2 aminoacids in a window which have a higher Blosum62 score than a given cutoff (here: 1)
 * 
 */
template <typename TSeq, typename TVertex, typename TVertexArray>
TGraph fastCreateVertices ( TSeq const & seq, unsigned const & winLen, unsigned const & maxIndelTandem, unsigned const & maxIndelCryptic, TVertex & vertex, TVertexArray & vertexArray, int const & originalImp ) {
    TGraph gbaGraph;
    int cutoff = 1;
    int counter = 0;
    //dummySource :
    int firstDumPos = 0;
    int secDumPos = 0;
    
    String<int> alphabetFreq = iniFreqVector2;

    for ( int i = 0; i < winLen; ++i ) {
        char currentLetter = seq[i];
        ++alphabetFreq[(int)(currentLetter)];
    }
    
    createVertexArray ( vertex, vertexArray, counter, firstDumPos, secDumPos);
    addVertex(gbaGraph);
    
    ++counter;
    
    for ( unsigned i = 0; i < winLen-1; ++i ) {
        char firstLetter = seq[i];
        
        for ( unsigned x = i+1; x < winLen; ++x ) {
            char secondLetter = seq[x];

            //std::cout << firstLetter << " " << secondLetter << " " <<  score ( blosum62Matrix, firstLetter, secondLetter ) << std::endl;
            if ( score ( blosum62Matrix, firstLetter, secondLetter ) > cutoff ) {
                
                if ( checkSameWindowByChance ( repeatMatrix, nonRepeatMatrix, firstLetter, secondLetter, alphabetFreq[firstLetter], winLen )) {
                    int firstPositionOfTuple = i+1;
                    int secondPositionOfTuple = x+1;
                    createVertexArray ( vertex, vertexArray, counter, firstPositionOfTuple, secondPositionOfTuple);
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

            //TODO in Kahveci's code, switched secondLetter and firstLetter order?
            char secondLetter = seq[i+x];
            char firstLetter = seq[i+winLen-1];
            if (originalImp==1) {    
                firstLetter = seq[i+x];
                secondLetter = seq[i+winLen-1]; 
            }
            

            if ( score ( blosum62Matrix, firstLetter, secondLetter ) > cutoff ) {
                //std::cout << firstLetter << " " << secondLetter << " " <<  score ( blosum62Matrix, firstLetter, secondLetter ) << std::endl;
            //TODO in Kahveci's code, switched secondLetter and firstLetter order?
                if ( checkSameWindowByChance ( repeatMatrix, nonRepeatMatrix, firstLetter, secondLetter, alphabetFreq[firstLetter], winLen )) {
                    int firstPositionOfTuple = i+x+1;
                    int secondPositionOfTuple = i+winLen;
                    
                    createVertexArray ( vertex, vertexArray, counter, firstPositionOfTuple, secondPositionOfTuple);
                    addVertex ( gbaGraph );

                    ++counter;
                } else {
                    //std::cout << "Symbols in window by Chance: " << i+x << " - " << i+winLen-1 << std::endl;
                }
            }
        }
    }
    
    //for (int i = 0; i < length(vertexArray); ++i) {
    //    printSeqAnString(vertexArray[i]);
    //}
    //std::cout << "Number of vertices in the graph: " <<numVertices(gbaGraph) << std::endl;
    
    return gbaGraph;
}


/* Calculates shannon entropy
 * 
 */
long double shanEntropy ( unsigned const & L, String<int> const & aaQuantity ) {
    long double K2KV = 0;
    
    for ( unsigned i = 0; i < 27; ++i ) {    
        long double conv=(long double)aaQuantity[i]/L;
        if ( aaQuantity[i] != 0 ) {
            K2KV += -log2(conv) * conv;
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

/*
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

/*
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
 * TODO Precalculate factorials and look for the results in a list b4 running the tool to save time rather than calculating them
 * thousands of times
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

/* reverse the order of array entries, e.g.: 0 0 0 2 4 => 4 2 0 0 0
 * 
 */
void reverse(int arr[], int count)
{
   int temp;
   for (int i = 0; i < count/2; ++i)
   {
      temp = arr[i];
      arr[i] = arr[count-i-1];
      arr[count-i-1] = temp;
   }
}

/* Calculates nat Log of colourings, the number of compositions in any complexity state
 * See equation 1 of Wootton and Federhen (March 1993)
 * 
 * */
double lnAss (String<int> const & freqVector) {
    long double ans;
    
    String<int> sortedVector = freqVector;
    
    int tmpl=length(freqVector);
    int* v_copy_tmp = new int[tmpl];
    std::copy(&sortedVector[0],&sortedVector[tmpl],&v_copy_tmp[0]);
    quickSort(v_copy_tmp, 0, tmpl-1);
    std::reverse<int*>(&v_copy_tmp[0],&v_copy_tmp[tmpl]);
    std::copy(&v_copy_tmp[0],&v_copy_tmp[tmpl],&sortedVector[0]);
    
    
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
    long double subSeq = ((long double) len) * (log((long double) alphSize));
    ans1 = lnAss (freqVector);
    
    if (ans1 > -100000.0) {
        ans2 = lnPerm (freqVector, len);
    } else {
        //TODO
        std::cout << "ERROR: INSERT ERRORLOG IN FUNCTION GETPROB" << std::endl;
    }
    ans = ans1 + ans2 - subSeq;
    
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
            String<AminoAcid> window = infix (tempSeq, winStartPos+count, winEndPos+count+1);
            
            if (length(window) <= 1) {
                break;
            }
            
            String<int> freqVector = iniFreqVector;
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

/*
 * 
 */
String<String<int> > getSourceVertexForSubGraph ( TGraph const & gbaGraph, int const & numComponents, String<unsigned> components )
{
    
    String<String<int> > sourceVertexList;
    resize (sourceVertexList, numComponents);
    int component;
    
    int vertexCount = numVertices ( gbaGraph );
    
    for (int i = 0; i < numComponents; ++i) {
        sourceVertexList[i] = 0;
        appendValue(sourceVertexList[i], 0); //
    }

    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it ( gbaGraph );
    while (!atEnd(it)) {
        if ( inDegree ( gbaGraph, getValue(it) ) == 0 ) {
            component = getProperty(components, getValue(it)); //this is the ID of the subgraph for this Vertex
            sourceVertexList[component][0] = component;
            appendValue(sourceVertexList[component], getValue(it));
        }
        goNext(it);
    }
    


    return sourceVertexList;
}

/* getPaths() finds all vertex IDs of a path longer than a cutoff 
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
        appendValue ( longPaths, source );
    } else if ( getProperty ( predecessor, v ) == getNil<typename VertexDescriptor<Graph<Directed<TCargo> > >::Type>() ) {
        std::cout << "No path from " << source << " to " << v << " exists.";

    } else {
        getPaths ( g,predecessor, source, getProperty ( predecessor, v ), longPaths );
        appendValue ( longPaths, v );
    }

}

/* getLongestPathsPerSubgraph() returns the longest path in a directed graph (reversed dijkstra algorithm, edgeValue * -1)
 * We use all vertices without any incoming edges as source vertices.
 * The original implementation (Kahveci, Li, 2006) extract every connected subgraph.
 * 
 */
String<int> getLongestPathsPerSubgraph(TGraph const & gbaGraph) {
    String<String<int> > vertexIdComponents;
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;

    TVertexIterator it ( gbaGraph );
    String<unsigned> components;
    unsigned numComponents = 0;
    String<int> predMap;
    String<int> distMap;
    String<int> weightMap;
    String<int> weights;
    String<int> longPaths;
    StringSet<String<int> > longPathsEachSubgraph;
    
    resize ( weights, numEdges(gbaGraph));
    
    int component = 0;
    
    for ( int i = 0; i < numEdges ( gbaGraph ); ++i ) {
        weights[i] = -1;
    }

    assignEdgeMap ( weightMap, gbaGraph, weights );
    
    numComponents = weaklyConnectedComponents(components, gbaGraph);

    String<String<int> >sourceVertices = getSourceVertexForSubGraph ( gbaGraph, numComponents, components );
    
    resize (vertexIdComponents, numComponents);
    for (int i = 0; i < numComponents; ++i) {
        vertexIdComponents[i] = sourceVertices[i][0];
        appendValue(vertexIdComponents[i], 0); //minimumDist
        appendValue(vertexIdComponents[i], 0); //node With longest Distance
        appendValue(vertexIdComponents[i], 0); // end node with longest Distance
        
    }
    resize (longPathsEachSubgraph, numComponents);
    
    //this is the most expensive part in GBA for big sequences
    //potential for optimization !!
    while (!atEnd(it)) {
        
        if ( outDegree ( gbaGraph, getValue(it) ) == 0 ) {
           
        
            component = getProperty(components, getValue(it)); //this is the ID of the subgraph for this Vertex
            
            for (int i = 2; i < length(sourceVertices[component]); ++i) {
                dagShortestPath ( predMap, distMap, gbaGraph, sourceVertices[component][i] ,weightMap );
                
                if ( getProperty ( distMap, getValue ( it ) ) < vertexIdComponents[component][1]) {
                    //std::cout << "path from: " << sourceVertices[component][i] << " to " << getValue(it) << " and distance is: " << getProperty(distMap, getValue(it));
                    //_printPath(gbaGraph, predMap, (TVertexDescriptor) sourceVertices[component][i], getValue(it));
                    vertexIdComponents[component][1] = getProperty(distMap, getValue(it)); //the longest distance
                    vertexIdComponents[component][2] = sourceVertices[component][i]; //saving the sourcenode from the longest path from component
                    vertexIdComponents[component][3] = getValue(it);
                    
                    getPaths ( gbaGraph, predMap, ( TVertexDescriptor ) vertexIdComponents[component][2], (TVertexDescriptor) vertexIdComponents[component][3], longPaths );    
                    
                    resize (longPathsEachSubgraph[component], length(longPaths));
                    longPathsEachSubgraph[component] = longPaths;
                    clear(longPaths);
                }
            }
        }
        goNext(it);
    }
    
    //for (int i = 0; i < length(sourceVertices); ++i) {
    //    printSeqAnString(sourceVertices[i]);
    //}
    
    String<int> concatLongPaths = concat(longPathsEachSubgraph);
    //std::cout << "long paths: " << std::endl;
    //printSeqAnString(concatLongPaths);
    return concatLongPaths;
}

/* datastructure for the sampleLenPercentage file (post-processing, gba)
 * 
 */

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

/* 
 * 
 */
String<String<double> > readSampleLenPercentage (String<char> const & pathToFile) {
    String<char> sampleLenPerFile = getAbsolutePath (toCString(pathToFile));
    
    String<String<double> > sampLenPer;
    String<double> grp;
    ifstream readFile(toCString(sampleLenPerFile));
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


int readSeqFromFastaFile ( String<char> const & pathToFile, StringSet<String<AminoAcid> > & seqSet, StringSet<String<char> > & seqIds)
{
    String<char> seqFileName = getAbsolutePath ( toCString ( pathToFile ) );

    SeqFileIn seqFileIn;
    if ( !open ( seqFileIn, toCString ( seqFileName ) ) ) {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

    try {
        readRecords ( seqIds, seqSet, seqFileIn );

    } catch ( Exception const & e ) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

void writeRegsIntoFile (String<String<String<int> > > const & segments, StringSet<String<char> > seqIds, String<char> const & pathToFile) {
    String<char> seqFileName = getAbsolutePath ( toCString ( pathToFile ) );
    
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
 
 /* SeqAn String AminoAcid contains 27 characters, because "wildcards" like B for Aspartic Acid or Asparagine is included.
  * in order to get the best calculations when computing entropies we need to extract those wildcards
  */
bool areRegularAminoAcids (String<AminoAcid> AA) {
    char firstA = AA[0];
    char secA = AA[1];    
    
    if (firstA == (char)'B' || firstA == (char)'J' || firstA == (char)'Z' || firstA == (char)'X' || firstA == (char)'O' || firstA == (char)'U' || firstA == (char)'*' ) {
        return false;
    }
    if (secA == (char)'B' || secA == (char)'J' || secA == (char)'Z' || secA == (char)'X' || secA == (char)'O' || secA == (char)'U' || secA == (char)'*') {
        return false;
    }
    return true;
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
            if (areRegularAminoAcids(lets1) && areRegularAminoAcids(lets2)) {
                twoLetScore = twoLetScore + score(blosum62Matrix, lets1[0], lets2[0]);
                twoLetScore = twoLetScore + score(blosum62Matrix, lets1[1], lets2[1]);
                array[i][j] = twoLetScore / 2;
                
            } else {
                array[i][j] = 0;
            }
            
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


template <typename tMat, size_t rows, size_t cols>
double computeTwoLetterEntropy(String<AminoAcid> const & subString, tMat (&twoLetterScoringMatrix)[rows][cols]) {
    
    double entropy = 0.;
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
                singleProb = singleProb + (twoLetterFrequencies[j] * twoLetterScoringMatrix[i][j]);
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
    
    return entropy;
}

/* Extend current blocks to left or right in single steps 
 * if the extended region has a lower entropy, we replace the old region
 */  
String<String<int> > extend (int startPos, int endPos, int limit, string const direction, String<AminoAcid> const & seq) {
    double comp1 = 0;
    //std::cout << "startPos: " << startPos << " endpos: " << endPos << " limit: " <<limit << std::endl;
    
    double comp2 = 0;
    int pointer = 0;
    int startDecPos = startPos;
    int endDecPos = endPos;
    
    String<AminoAcid> extReg = infix(seq, startPos-1, endPos);
    String<String<int> > decRegs;
    clear (decRegs);
    
    
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
                    endDecPos = pointer+2;
                }    
            } else if (dec) {
                dec = false;
                startDecPos = pointer + 2;
                if (comp1 < comCut) {
                    String<String<int> > found = startDecPos;
                    appendValue(found[0], endDecPos);
                    append(found, decRegs);
                    resize (decRegs, length(found));
                    decRegs = found;
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
                        String<String<int> > found = startDecPos;
                        appendValue(found[0], endDecPos);
                        append(found, decRegs);
                        resize (decRegs, length(found));
                        decRegs = found;
                        
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
                
                String<String<int> > found = startDecPos;
                appendValue(found[0], endDecPos);
                append(found, decRegs);
                resize (decRegs, length(found));
                decRegs = found;
                
            }
        }
        if (length(decRegs) == 0) {
            
            //String<int> addToDecRegs = startPos;
            //appendValue(addToDecRegs, endPos);
            //appendValue(decRegs, addToDecRegs);
            
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
                    startDecPos = pointer-1;
                }    
            } else if (dec) {
                dec = false;
                endDecPos = pointer -1;
                if (comp1 < comCut) {
                    String<String<int> > found = startDecPos;
                    appendValue(found[0], endDecPos);
                    append(decRegs, found);
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
                        
                        String<String<int> > found = startDecPos;
                        appendValue(found[0], endDecPos);
                        append(decRegs, found);                      
                    }
                }
                ++pointer;
            }
            
            if ((pointer == limit) && (dec)) {
                endDecPos = limit - 1;
                
                String<String<int> > found = startDecPos;
                appendValue(found[0], endDecPos);
                append(decRegs, found);
   
            }
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

/* returns regions blockwise by position in the sequence from a list e.g. 4-7, 9-10 from 4 5 6 7 9 10 
 * 
 */
String<String<int> > genBlocksFromList(String<int> const & blockList) {
    String<int> tmpBlockList = blockList;
    String<String<int> > lcrBlocks;
    int pos = 0;
    int counter = 0;
    int startPosLcrBlock = tmpBlockList[pos];    
    
    for (unsigned i = 0; i < length(tmpBlockList) ; ++i) {
        if (tmpBlockList[i]+1 !=tmpBlockList[i+1] ) {
            int endPosLcrBlock = tmpBlockList[i];
            if (endPosLcrBlock - startPosLcrBlock > 1) {    
                appendValue(lcrBlocks, startPosLcrBlock);
                appendValue(lcrBlocks[counter], endPosLcrBlock);
                ++counter;
            }
            //std::cout << startPosLcrBlock << "-" << endPosLcrBlock << std::endl;
            startPosLcrBlock = tmpBlockList[pos+1];
            
        }
        ++pos;
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

/*
 * 
 */
bool shareLetter (String<AminoAcid> const & str1, String<AminoAcid> const & str2) {
    bool shared = false;
    
    for (int i = 0; i < length(str1); ++i) {
        for (int j = 0; j < length(str2); ++j) {
            if (str1[i] == str2[j]) {
                shared = true;
                return shared;
            }
        }
    }
    
    return shared;
}


/*
 * 
 */
bool checkContribution (String<int> const & currentBlock, String<String<int> > const & decRegs, String<AminoAcid> const & seq) {
    bool contributed = false;
    String<String<int> > regs = decRegs;
    String<int> block;
    int i = 0;
    int len = length(regs);
    while ((i < len) && (!contributed)) {
        block = regs[i];
        int start = block[0];
        int end = block[1];
        String<AminoAcid> subBlock = infix(seq, start, end+1);
        String<AminoAcid> curBlock = infix(seq, currentBlock[0], currentBlock[1]+1);
        contributed = shareLetter(curBlock, subBlock);
        ++i;
    }
    return contributed;
}


/* working on LCR blocks (extending while considering limits to the left and right)
 * limits are the end of the sequence or start/end positions of found blocks
 */
String<String<int> > processBlocks (String<String<int> > const & blocks, String<AminoAcid> const & seq) {
    String<String<int> > tmpBlocks = blocks;
    String<String<int> > frontLcrs;
    String<String<int> > backLcrs;
    String<String<int> > lcrs;
    String<int> currentBlock;
    resize (currentBlock, 2);
    bool isFirstBlock = true;
    String<int> saveExtendedBlock;

    String<int> tmpBlock;
    int limit = 0;
    int startPos, endPos = 0;
    
    int counter = 0;
    while (length(tmpBlocks) != 0) {
        int lcrBlockStart = 0;
        int lcrBlockEnd = 0;    
        clear(frontLcrs);
        clear(backLcrs);
        bool extendToLeft = true;
        bool find = false;
        int lastElement = length(lcrs)-1;
        if (length(lcrs)!=0) {
            
            tmpBlock = lcrs[lastElement];
            
            lcrBlockEnd = tmpBlock[1];
            while ((!find)&& (length(tmpBlocks)!=0)) {

                currentBlock = tmpBlocks[0];
                
                erase(tmpBlocks, 0);
                
                startPos = currentBlock[0];
                endPos = currentBlock[1];
                
                if (startPos < lcrBlockEnd) {
                    if (endPos > lcrBlockEnd) {
                        if ((endPos-lcrBlockEnd) >= 3) {
                            
                            startPos = lcrBlockEnd + 1;
                            extendToLeft = false;
                            find = true;
                        }
                    }
                } else {
                    find = true;
                }
            }
        } else {
            
            currentBlock = tmpBlocks[0];
            
            erase(tmpBlocks, 0);
            
            startPos = currentBlock[0];
            endPos = currentBlock[1];
            find = true;
        }
        
        
        if (find==true){
            if (isFirstBlock==true) {
                limit = -1;
                isFirstBlock=false;
                frontLcrs = extend (startPos, endPos, limit, "left", seq);
            } else if (extendToLeft){
                
                limit = lcrBlockEnd -1;
                frontLcrs= extend (startPos, endPos, limit, "left", seq);
            }
            
            limit = length(seq)+1;
            backLcrs = extend(startPos, endPos, limit, "right", seq);
            
            double com = 0.;
            int cbStart = currentBlock[0];
            int cbEnd = currentBlock[1];
            String<AminoAcid> tmpSeq = infix(seq, cbStart, cbEnd+1);
            com = computeTwoLetterEntropy(tmpSeq, twoLetterScoringMatrix);
            bool contributed = false;
            
            if (length(frontLcrs) != 0) {
                tmpBlock = frontLcrs[0];
                lcrBlockStart = tmpBlock[0];
                if (com > comCut) {
                    contributed = checkContribution(currentBlock, frontLcrs, seq);
                    if (!contributed) {
                        lcrBlockEnd = startPos - 1;
                        
                    } else {
                        lcrBlockEnd = endPos;
                    }
                } else {
                    lcrBlockEnd = endPos;
                }
                
                String<int> addToLcrs = lcrBlockStart;
                appendValue(addToLcrs, lcrBlockEnd);
                appendValue (lcrs, addToLcrs);
            }
            
            bool combine = false;
            
            if ((!contributed) && (com > comCut)) {
                contributed = checkContribution(currentBlock, backLcrs, seq);
            
                if (!contributed) {
                        lcrBlockStart = endPos + 1;
                } else {
                    if (length(frontLcrs) != 0) {
                        combine = true;
                    }
                }
            } else if (length(frontLcrs) != 0) {
                combine = true;
            }
            
            if (length(backLcrs) != 0) {
                int lastElemBackLcr = length(backLcrs)-1;
                tmpBlock = backLcrs[lastElemBackLcr];
                lcrBlockEnd = tmpBlock[1];
                
                if (combine) {
                    limit = length(lcrs);
                    tmpBlock = lcrs[limit-1];
                    
                    erase(lcrs, limit-1);
                    lcrBlockStart = tmpBlock[0];
                    
                } else {
                    if (com < comCut ) {
                        lcrBlockStart = startPos;
                    } else if (contributed) {
                        lcrBlockStart = startPos;
                    } else {
                        lcrBlockStart = endPos + 1;
                    }
                }
                
                String<int> addToLcrs = lcrBlockStart;
                appendValue(addToLcrs, lcrBlockEnd);
                appendValue (lcrs, addToLcrs);
                
                
            } else {
                if ((length(frontLcrs) == 0) && (!contributed) && (com < comCut)) {
                    
                    String<int> addToLcrs = currentBlock[0];
                    appendValue(addToLcrs, currentBlock[1]);
                    appendValue (lcrs, addToLcrs);
                
                }
            }
        }
    }
    
    return lcrs;
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
        
        if ((curBlockStart <= nextBlockStart) && (nextBlockStart-1 <= curBlockEnd)) {
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
        com = com/(length(str)-1); //normalize
        
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
        com = com/(length(str)-1); //normalize
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
        
        auto src1Pos = toSourcePosition(row(align,0),0);
        auto src1End = toSourcePosition(row(align,0),len-1);
        auto src2Pos = toSourcePosition(row(align,1),0);
        auto src2End = toSourcePosition(row(align,1),len-1);
        
        int similarity = 0;
        
        
        for (unsigned i = 0; i < len; ++i) {
            AminoAcid firstChar = row(align,0)[i];
            AminoAcid secChar = row(align,1)[i];
            
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
 * mark == -1 means front, mark == -2 means back
 */
StringSet<String<int> > checkAdjBlock (int const & start1, int const & end1, String<int> adjBlock, String<AminoAcid> const & seq, double const & cCut, int const & mark) {    
    double com = 0;
    StringSet<String<int> > result;
    String<AminoAcid> seq1 = infix(seq, start1-1, end1);
    int start2 = adjBlock[0];
    int end2 = adjBlock[1];
    //TODO better?  infix(seq, start2-2, ...);
    String<AminoAcid> seq2 = infix(seq, start2-1, end2);
    
    StringSet<String<int> > aliPos;
    if (mark == -1) {
        aliPos = findAlignment (seq2, seq1);
    } else {
        aliPos = findAlignment (seq1, seq2);
    }
    
    /* 
    std::cout << "AliPos: " << std::endl;
    for (int i = 0; i < length(aliPos); ++i) {
        printSeqAnString(aliPos[i]);
    }
    */
    
    
    if (length(aliPos) != 0) {
        int aliStart2 = aliPos[0][0];
        int aliEnd2 = aliPos[0][1];
            
        int aliStart1 = aliPos[1][0];
        int aliEnd1 = aliPos[1][1];
        
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

            
            addToResult (result, left);
            
            append1 = start2 + aliStart1 ;
            append2 = start2 + aliEnd1 ;
            appendValue (resultPrep, append1);
            appendValue (resultPrep, append2);
            appendValue (result, resultPrep);
            clear (resultPrep);
            
            left = checkLeftRegs (aliStart1, aliEnd1, start2, end2, seq, cCut);
            
            addToResult (result, left);
            
        }
        
    }
    
    return result;
}


/* checks if there are deletable characters in a region 
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
    
    //std::cout << "--------------------------FILTERING PROCESS----------------------------"<<std::endl;
    int seqLen = length(seq);
    int len = length(mergedBlocks);
    String<double> complexityList ;
    String<String<int> > tmpBlocks = mergedBlocks;
    String<String<double> > sampLenPer = readSampleLenPercentage("../LCR-finder/matrices/sampledLenRepPer");
    double comp = 0.;
    double max = -9999999;
    int it = 0;
    
    if (len != 1) {
        for (unsigned i = 0; i < len; ++i) {
            int start = tmpBlocks[i][0];
            int end = tmpBlocks[i][1];
            String<AminoAcid> str = infix(seq, start-1, end);
            
            comp = computeTwoLetterEntropy(str, twoLetterScoringMatrix);
            comp = comp/(length(str)-1); //normalize
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
            per = sampLenPer[lineNr][3];
            
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
            
            it = 0;
            max = -9999999;
            while (it < compLen) {
                double comp = complexityList[it];
                if (comp > max) {
                    max = comp;
                    maxIndex = it;
                }
                ++it;
            }
            
            cCut = complexityList[maxIndex];
            erase(complexityList, maxIndex);
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
            comp = comp/(length(str)-1); //normalize
            
            if (comp >= cCut ) {
                result = checkDeletability (tmpBlocks, it, seq, cCut);
                
                int rSize = length(result);
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
                    //std::cout << "removed (high complexity)" << std::endl;
                    erase (tmpBlocks, it);
                    ++lineNr;
                    len = length(tmpBlocks);
                }
            } else {
                ++it;
            }
        }
    }
    
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
        
        String<int> freqVector = iniFreqVector;
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



/* SEG algorithm starts here
 * 
 */
String<String<int> > seg ( String<AminoAcid> const & seq, unsigned const & winLen, float const & K2A, float const & K2B )
{
    
    if ( winLen > length ( seq ) | K2A > K2B ) {
        std::cerr << "Input-parameters are not allowed." << std::endl;
        return 0;
    }
    
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
        std::cout << std::endl;
    }
    
    return mergedRegions;
}


/* GBA algorithm starts here
 * 
 */
String<String<int> > gba (String<AminoAcid> const & seq, int const winLen, int const maxIndelTandem, int const maxIndelCryptic, int const & originalImp) {
    
    String<int> vertex;
    String<String<int> > vertexArray;

    if ( winLen > length ( seq ) | maxIndelTandem >= winLen | maxIndelCryptic >= winLen ) {
        std::cerr << "Input-parameters are not allowed." << std::endl;
        return 0;
    }
    
    TGraph gbaGraph = fastCreateVertices ( seq, winLen, maxIndelTandem, maxIndelCryptic, vertex, vertexArray, originalImp );
    gbaGraph = createEdges ( gbaGraph, vertexArray, seq );
    String<int> longPathVertices = getLongestPathsPerSubgraph(gbaGraph);
    
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
    
    return filteredRegions;
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
    int originalImp;
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
    addOption(parser, seqan::ArgParseOption("oi", "originalImp", "Choose whether you want to use the original Implementation of creating Vertices or the modified version.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "originalImp", 1);
    
    
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
    getOptionValue(options.originalImp, parser, "originalImp");
    
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
    int  originalImp = options.originalImp;
    
    
    
    
      
    if (alg == "gba" or alg == "GBA") {
        if (readSeqFromFastaFile ( fileIn, seqSet, seqIds )) {
            return 1;
        }
        std::cout << "------------------------------- executing GBA... -----------------------------" << std::endl;
        computeTwoLetterScoringMatrix(twoLetterScoringMatrix, twoGramAlphabet);

        normalizeScoringMatrix(twoLetterScoringMatrix);
        SEQAN_OMP_PRAGMA (parallel for)
        for (unsigned i = 0; i < length(seqSet); ++i) { 
            String<AminoAcid> seq = seqSet[i];
            segs = gba(seq, windowLength, maxIndelTandem, maxIndelCryptic, originalImp);
            appendValue (segments, segs);
            
        }
        
    } else if (alg == "seg" or alg == "SEG") {
        if (readSeqFromFastaFile ( fileIn, seqSet, seqIds )) {
            return 1;
        }
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
      
    return 0;
}
