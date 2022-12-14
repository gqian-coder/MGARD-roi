#ifndef MGARD_X_MESSGAE_HH
#define MGARD_X_MESSGAE_HH

#include <iostream>
#include <sstream>
#include <string>

using std::string;

namespace mgard_x {
namespace log {

extern const string log_null;
extern const string log_err;
extern const string log_dbg;
extern const string log_info;
extern const string log_warn;
extern const string log_time;

// https://stackoverflow.com/a/26080768/8740097
template <typename T> void build(std::ostream &o, T t);

template <typename T, typename... Args>
void build(std::ostream &o, T t, Args... args);

template <typename... Args> void print(string log_head, Args... args);

} // namespace log

} // namespace mgard_x

#endif // FORMAT_HH
