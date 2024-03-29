/*! Definition of primitives for a text based input/output.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#include "TauTriggerTools/Common/interface/TextIO.h"

#include <cmath>
#include <iomanip>
#include <unordered_set>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp>
#include "TauTriggerTools/Common/interface/exception.h"

namespace analysis {

std::string RemoveFileExtension(const std::string& file_name)
{
    return boost::filesystem::change_extension(file_name, "").string();
}

std::string GetFileNameWithoutPath(const std::string& file_name)
{
    const size_t lastindex = file_name.find_last_of("/");
    if(lastindex == std::string::npos)
        return file_name;
    else
        return file_name.substr(lastindex+1);
}

std::vector<std::string> SplitValueList(const std::string& _values_str, bool allow_duplicates,
                                        const std::string& separators, bool enable_token_compress)
{
    std::string values_str = _values_str;
    std::vector<std::string> result;
    if(enable_token_compress)
        boost::trim_if(values_str, boost::is_any_of(separators));
    if(!values_str.size()) return result;
    const auto token_compress = enable_token_compress ? boost::algorithm::token_compress_on
                                                      : boost::algorithm::token_compress_off;
    boost::split(result, values_str, boost::is_any_of(separators), token_compress);
    if(!allow_duplicates) {
        std::unordered_set<std::string> set_result;
        for(const std::string& value : result) {
            if(set_result.count(value))
                throw exception("Value '%1%' listed more than once in the value list '%2%'.") % value % values_str;
            set_result.insert(value);
        }
    }
    return result;
}

std::vector<std::string> ReadValueList(std::istream& stream, size_t number_of_items, bool allow_duplicates,
                                       const std::string& separators, bool enable_token_compress)
{
    const auto stream_exceptions = stream.exceptions();
    stream.exceptions(std::istream::goodbit);
    try {
        std::vector<std::string> result;
        std::unordered_set<std::string> set_result;
        const auto predicate = boost::is_any_of(separators);
        size_t n = 0;
        for(; n < number_of_items; ++n) {
            std::string value;
            while(true) {
                const auto c = stream.get();
                if(!stream.good()) {
                    if(stream.eof()) break;
                    throw exception("Failed to read values from stream.");
                }
                if(predicate(c)) {
                    if(!value.size() && enable_token_compress) continue;
                    break;
                }
                value.push_back(static_cast<char>(c));
            }
            if(!allow_duplicates && set_result.count(value))
                throw exception("Value '%1%' listed more than once in the input stream.") % value;
            result.push_back(value);
            set_result.insert(value);
        }
        if(n != number_of_items)
            throw exception("Expected %1% items, while read only %2%.") % number_of_items % n;

        stream.clear();
        stream.exceptions(stream_exceptions);
        return result;
    } catch(exception&) {
        stream.clear();
        stream.exceptions(stream_exceptions);
        throw;
    }
}

StVariable::StVariable() : value(0), error_up(0), error_low(0) {}
StVariable::StVariable(double _value, double _error_up, double _error_low) :
    value(_value), error_up(_error_up), error_low(_error_low) {}

int StVariable::precision_up() const
{
    return error_up != 0.
            ? static_cast<int>(std::floor(std::log10(error_up)) - number_of_significant_digits_in_error + 1)
            : max_precision;
}

int StVariable::precision_low() const
{
    return error_low != 0.
            ? static_cast<int>(std::floor(std::log10(error_low)) - number_of_significant_digits_in_error + 1)
            : max_precision;
}

int StVariable::precision() const { return std::max(precision_up(), precision_low()); }

int StVariable::decimals_to_print_low() const { return std::max(0, -precision_low()); }
int StVariable::decimals_to_print_up() const { return std::max(0, -precision_up()); }
int StVariable::decimals_to_print() const { return std::min(decimals_to_print_low(), decimals_to_print_up()); }

std::string StVariable::ToLatexString() const
{
    const ValueType ten_pow_p = std::pow(10.0, precision());
    const ValueType value_rounded = std::round(value / ten_pow_p) * ten_pow_p;
    const ValueType error_up_rounded = std::ceil(error_up / ten_pow_p) * ten_pow_p;
    const ValueType error_low_rounded = std::ceil(error_low / ten_pow_p) * ten_pow_p;

    std::ostringstream ss;
    ss << std::setprecision(decimals_to_print()) << std::fixed;
    if(error_up == 0 && error_low == 0)
        ss << value_rounded<< "^{+0}_{-0}";
    else if(!std::isnan(error_low))
        ss << value_rounded<< "^{+" << error_up_rounded << "}_{-" << error_low_rounded << "}";
    else if(std::isnan(error_low)) {
        ss << value_rounded << " \\pm ";
        if(error_up == 0)
            ss << "0";
        else
            ss << error_up_rounded;
   }

   return ss.str();
}

} // namespace analysis
