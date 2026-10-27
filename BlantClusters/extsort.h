#pragma once
#include <cstddef>
#include <memory>
#include <string>

// The numeric value of the record's first whitespace-separated field.
double leadingNumber(const std::string& record);

// Sorts newline-free text records that begin with a whitespace-separated
// double field, descending, ties broken by the whole record ascending.
// Spills runs of at most bufBytes to runPrefix + ".runN", then merges them.
// bufBytes counts each record's text plus its std::string header; allocator
// slack can still put real use somewhat above it.
class ExternalSorter {
public:
    ExternalSorter(std::string runPrefix, size_t bufBytes);
    ~ExternalSorter();
    ExternalSorter(const ExternalSorter&) = delete;
    ExternalSorter& operator=(const ExternalSorter&) = delete;

    void add(std::string record);  // record must contain no '\n'
    long size() const;             // records added so far
    void finish();                 // flush the last run, start merging
    bool next(std::string& out);   // valid after finish(); false at end
    void discard();                // close and unlink the run files

private:
    struct Impl;
    std::unique_ptr<Impl> p;
};
