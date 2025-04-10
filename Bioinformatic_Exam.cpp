#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>
#include <numeric>
using namespace std;

bool hasEnding(string const& fullString, string const& ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    }
    else {
        return false;
    }
}
bool IsInBounds(const int& value, const int& low, const int& high) {
    return !(value < low) && !(high < value);
}

// Function to read sequences from a FASTA file
vector<string> readFASTAFile(string& filename) {
    ifstream file(filename);
    while (!file.is_open()) {
        cerr << "Error: Could not open " << filename << ". Please enter a valid filename!" << endl;
        cin >> filename;
        if (!hasEnding(filename, ".fasta")) {
            filename += ".fasta";
            cout << "Note: .fasta extension was added. Using: " << filename << endl;
        }
        file = ifstream(filename);
    }

    vector<string> sequences;
    string line;
    string seq;

    while (getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!seq.empty()) {
                sequences.push_back(seq);
                seq = "";
            }
        }
        else {
            seq += line;
        }
    }

    // Incase the FASTA-file does not contain two sequences
    if (!seq.empty()) sequences.push_back(seq);
    if (sequences.size() < 2) {
        cerr << "Error: File did not contain at least 2 sequences!" << endl;
        exit(1);
    }
    return sequences;
}

// Ensure that scoring parameters only takes numbers
int getValidInteger(const string& prompt) {
    double input;
    double test = 0;
    int value;
    while (true) {
        cout << prompt;
        cin >> input;
        if (cin.fail()) {
            cin.clear(); // Clear the error flag
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); // Discard invalid input
            cout << "Invalid input. Please enter a number." << endl;
        }
        else {
            if (modf(input, &test) != 0) {
                value = (int)round(input); //round of the input into the nearest integer value
                cout << "input was not a whole number. It will be rounded to " << value << endl;
            }
            else {
                value = (int)input;
            }
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); // Discard extra input (e.g., '5abc')
            return value;
        }
    }
}
//ensure that the input is in the permitted values
int getValidIntegerWithinPermittedValues(const string& prompt, vector<int> permittedValues) {
    double input;
    double test = 0;
    int value;
    while (true) {
        cout << prompt;
        cin >> input;
        if (cin.fail()) {
            cin.clear(); // Clear the error flag
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); // Discard invalid input
            cout << "Invalid input. Please enter a number." << endl;
        }
        else {
            if (modf(input, &test) != 0) {
                value = (int)round(input);
                cout << "input was not a whole number. It will be rounded to " << value << endl;
            }
            else {
                value = (int)input;
            }
            if (std::find(permittedValues.begin(), permittedValues.end(), value) == permittedValues.end()) {
                cerr << "Error: invalid input, please enter a permissiable number" << endl;
                return getValidIntegerWithinPermittedValues(prompt, permittedValues);
            }
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); // Discard extra input (e.g., '5abc')
            return value;
        }
    }
}

// Print the score matrix calculations
void printScoreMatrix(const vector<vector<int>>& dp,
    const string& seq1, const string& seq2,
    ofstream& outputFile) {
    int n = seq1.length();
    int m = seq2.length();

    // Determine the maximum width needed for any cell
    int max_width = 2; // Minimum width
    for (const auto& row : dp) {
        for (int val : row) {
            int width = to_string(val).length();
            max_width = max(max_width, width);
        }
    }
    int col_width = max_width + 1; // Each column width

    // Print the horizontal sequence header
    outputFile << string(col_width + 1, ' '); // Offset for first column
    for (int j = 0; j < m; j++) {
        outputFile << setw(col_width) << seq2[j];
    }
    outputFile << endl << endl;

    // Print the matrix with perfect alignment
    for (int i = 0; i <= n; i++) {
        // Print row label (vertical sequence)
        if (i > 0) {
            outputFile << seq1[i - 1] << " ";
        }
        else {
            outputFile << "  "; // Space for first row
        }

        // Print all values in the row
        for (int j = 0; j <= m; j++) {
            outputFile << setw(col_width) << dp[i][j];
        }
        outputFile << endl;
    }
}

// Function to perform Needleman-Wunsch Global Alignment
void needlemanWunsch(string seq1, string seq2, double match, double mismatch, double gap, ofstream& outputFile) {
    int n = seq1.length();
    int m = seq2.length();

    // Create DP matrix and traceback matrix
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));
    vector<vector<char>> trace(n + 1, vector<char>(m + 1, ' '));

    // Initialize the matrix
    for (int i = 0; i <= n; i++) {
        dp[i][0] = i * gap;
        trace[i][0] = 'U'; // Up
    }
    for (int j = 0; j <= m; j++) {
        dp[0][j] = j * gap;
        trace[0][j] = 'L'; // Left
    }

    // Fill the matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int scoreDiag = dp[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch);
            int scoreUp = dp[i - 1][j] + gap;
            int scoreLeft = dp[i][j - 1] + gap;
            dp[i][j] = max({ scoreDiag, scoreUp, scoreLeft });

            // Set traceback direction
            if (dp[i][j] == scoreDiag) {
                trace[i][j] = 'D'; // Diagonal
            }
            else if (dp[i][j] == scoreUp) {
                trace[i][j] = 'U'; // Up
            }
            else {
                trace[i][j] = 'L'; // Left
            }
        }
    }

    // Traceback to find the optimal alignment
    string align1, align2;
    int i = n, j = m;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && trace[i][j] == 'D') {
            align1 = seq1[i - 1] + align1;
            align2 = seq2[j - 1] + align2;
            i--;
            j--;
        }
        else if (i > 0 && trace[i][j] == 'U') {
            align1 = seq1[i - 1] + align1;
            align2 = '-' + align2;
            i--;
        }
        else {
            align1 = '-' + align1;
            align2 = seq2[j - 1] + align2;
            j--;
        }
    }

    // Output the alignment and score to the file
    outputFile << "========================================" << endl;
    outputFile << "     NEEDLEMAN-WUNSCH ALGORITHM" << endl;
    outputFile << "========================================" << endl;
    outputFile << "Sequence 1: " << seq1 << endl;
    outputFile << "Sequence 2: " << seq2 << "\n" << endl;
    outputFile << "Match Score: " << match << endl;
    outputFile << "Mismatch Score: " << mismatch << endl;
    outputFile << "Gap Penalty: " << gap << "\n" << endl;
    outputFile << "Alignment Score: " << dp[n][m] << endl;
    outputFile << "-----------------------------------------" << endl;
    outputFile << "Optimal Alignment: " << endl;
    outputFile << align1 << endl;
    outputFile << align2 << "\n" << endl;
    outputFile << "Scoring Matrix: " << endl;
    printScoreMatrix(dp, seq1, seq2, outputFile);
    outputFile << endl;
}

// Function to perform Smith-Waterman Local Alignment
void smithWaterman(string seq1, string seq2, int match, int mismatch, int gap, ofstream& outputFile) {
    int n = seq1.length();
    int m = seq2.length();

    // Create DP matrix and traceback matrix
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));
    vector<vector<char>> trace(n + 1, vector<char>(m + 1, ' '));

    int maxScore = 0;
    int maxI = 0, maxJ = 0;

    // Fill the matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int scoreDiag = dp[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch);
            int scoreUp = dp[i - 1][j] + gap;
            int scoreLeft = dp[i][j - 1] + gap;
            dp[i][j] = max({ 0, scoreDiag, scoreUp, scoreLeft });

            // Set traceback direction
            if (dp[i][j] == scoreDiag) {
                trace[i][j] = 'D'; // Diagonal
            }
            else if (dp[i][j] == scoreUp) {
                trace[i][j] = 'U'; // Up
            }
            else if (dp[i][j] == scoreLeft) {
                trace[i][j] = 'L'; // Left
            }

            // Track maximum score position
            if (dp[i][j] > maxScore) {
                maxScore = dp[i][j];
                maxI = i;
                maxJ = j;
            }
        }
    }

    // Traceback from the maximum score position
    string align1, align2;
    int i = maxI, j = maxJ;
    while (i > 0 && j > 0 && dp[i][j] != 0) {
        if (trace[i][j] == 'D') {
            align1 = seq1[i - 1] + align1;
            align2 = seq2[j - 1] + align2;
            i--;
            j--;
        }
        else if (trace[i][j] == 'U') {
            align1 = seq1[i - 1] + align1;
            align2 = '-' + align2;
            i--;
        }
        else {
            align1 = '-' + align1;
            align2 = seq2[j - 1] + align2;
            j--;
        }
    }

    // Output the alignment and score to the file
    outputFile << "========================================" << endl;
    outputFile << "     SMITH-WATERMAN ALIGNMENT" << endl;
    outputFile << "========================================" << endl;
    outputFile << "Sequence 1: " << seq1 << endl;
    outputFile << "Sequence 2: " << seq2 << "\n" << endl;
    outputFile << "Match Score: " << match << endl;
    outputFile << "Mismatch Score: " << mismatch << endl;
    outputFile << "Gap Penalty: " << gap << "\n" << endl;
    outputFile << "Local Alignment Score: " << maxScore << endl;
    outputFile << "-----------------------------------------" << endl;
    outputFile << "Optimal Local Alignment: " << endl;
    outputFile << align1 << endl;
    outputFile << align2 << "\n" << endl;
    outputFile << "Scoring Matrix:" << endl;
    printScoreMatrix(dp, seq1, seq2, outputFile);
    outputFile << endl;
}


int main() {
    string filename;
    string outputFilename;
    int match, mismatch, gap, choice;


    // Read sequences from the FASTA file and adds .fasta extension 
    cout << "Enter the filename with multiple sequence pairs (FASTA format): ";
    cin >> filename;
    if (!hasEnding(filename, ".fasta")) {
        filename += ".fasta";
        cout << "Note: .fasta extension was added. Using: " << filename << endl;
    }
    vector<string> sequences = readFASTAFile(filename);

    // Open the output file and adds .txt extension 
    cout << "Enter a name for the output file: ";
    cin >> outputFilename;
    if (!hasEnding(outputFilename, ".txt")) {
        outputFilename += ".txt";
        cout << "Note: .txt extension was added. Using " << outputFilename << endl;
    }

    ofstream outputFile(outputFilename);
    if (!outputFile.is_open()) {
        cerr << "Error: Could not open the output file!" << endl;
        return 1;
    }

    // Scoring parameters intake from user
    match = getValidInteger("Enter match score: ");
    mismatch = getValidInteger("Enter mismatch score: ");
    gap = getValidInteger("Enter gap penalty: ");
    choice = getValidIntegerWithinPermittedValues("Select alignment method (1 for Needleman-Wunsch, 2 for Smith-Waterman): ", { 1,2 });

    int seqLen = sequences.size();
    int firstSeq = 0;
    int secSeq = 1;

    vector<int> permittedValues(seqLen);
    // fills the vector from 1 to seqLen
    iota(permittedValues.begin(), permittedValues.end(), 1);

    // Incase there are more than two sequences in the FASTA-file
    if (seqLen > 2) {
        cout << "There are " << seqLen << " sequences in the file:\n";
        for (int i = 0; i < seqLen; i++) {
            cout << "Sequence " << (i + 1) << ": " << sequences[i].substr(0, 20) << "...\n";
        }
        std::ostringstream firstSeqString;
        firstSeqString << "Enter the number of the first sequence (1-" << seqLen << "): ";
        std::string prompt = firstSeqString.str();
        firstSeq = getValidIntegerWithinPermittedValues(prompt, permittedValues);

        std::ostringstream secSeqString;
        secSeqString << "Enter the number of the second sequence (1-" << seqLen << "): ";
        prompt = secSeqString.str();
        secSeq = getValidIntegerWithinPermittedValues(prompt, permittedValues);

        // Convert from 1-based to 0-based indexing
        firstSeq--;
        secSeq--;
    }

    if (choice == 1) {
        needlemanWunsch(sequences[firstSeq], sequences[secSeq], match, mismatch, gap, outputFile);
    }
    else if (choice == 2) {
        smithWaterman(sequences[firstSeq], sequences[secSeq], match, mismatch, gap, outputFile);
    }

    // Close the output file
    outputFile.close();
    cout << "Alignments and scores have been written to '" << outputFilename << "'." << endl;

    return 0;
}