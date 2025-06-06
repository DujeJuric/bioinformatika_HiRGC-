# bioinformatika_HiRGC

Project based on creating an own implementation of high performance referential genome compression algorithm.

### Credits

https://github.com/yuansliu/HiRGC

### Installation

g++ compiler is required.

From the project directory compile the source code:

```bash
g++ -std=c++17 -o hirgc_compr .\hirgc_compress.cpp -O2 -Wall
g++ -std=c++17 -o hirgc_decompr .\hirgc_decompress.cpp -O2 -Wall
```

### Data sets

2 data sets are provided in TestData folder, each containing refrence and target FASTA files

### Compression

Compressing a FASTA file:

```bash
./hirgc_compr <refrence FASTA file> <target FASTA file>
```

The result compressed txt file should be located in target folder from used data folder.

### Decompression

Decompressing a compressed txt file:

```bash
./hirgc_decompr <refrence FASTA file> <compressed txt file>
```

The result decompressed fasta file should be located in target folder from used data folder.
