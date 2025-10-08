#include "ExchangeWriter.h" 

void ExchangeWriter::write(const unordered_set<shared_ptr<Exchange>>& exchanges) {
	if (exchanges.empty()) return;
	for (const auto& e : exchanges) e->print(outFile);
}