#pragma once
#include <nlohmann/json.hpp>
#include <string>
#include <algorithm>
#include <type_traits>
#include <cctype>
#include <unordered_set>

namespace mrchem::json_utils {

using json = nlohmann::json;

// Strict-but-tolerant conversion to bool.
// Accepts true/false, 0/1, and strings like "on/off", "yes/no".
inline bool to_bool(const nlohmann::json& j, bool def = false) {
    if (j.is_boolean()) return j.get<bool>();
    if (j.is_number_integer()) return j.get<int>() != 0;
    if (j.is_number_float())   return j.get<double>() != 0.0;
    if (j.is_string()) {
        std::string s = j.get<std::string>();
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
        if (s == "true"  || s == "1" || s == "yes" || s == "on")  return true;
        if (s == "false" || s == "0" || s == "no"  || s == "off") return false;
        return def;
    }
    return def;
}

// Safely fetch a value of type T. For bool, uses to_bool.
// For numbers/strings, only returns if JSON type matches; else returns default.
template <class T>
inline T value_loose(const nlohmann::json& jparent, const char* key, const T& def) {
    if (!jparent.contains(key)) return def;
    const auto& j = jparent.at(key);
    if constexpr (std::is_same_v<T,bool>) {
        return to_bool(j, def);
    } else if constexpr (std::is_integral_v<T>) {
        if (j.is_number_integer()) return j.get<T>();
    } else if constexpr (std::is_floating_point_v<T>) {
        if (j.is_number()) return j.get<T>();
    } else if constexpr (std::is_same_v<T,std::string>) {
        if (j.is_string()) return j.get<std::string>();
    } else {
        try { return j.get<T>(); } catch (...) { /* fallthrough */ }
    }
    return def;
}

// ---------- New: input sanitizer to normalize booleans in the whole tree ----------

inline bool looks_like_bool_key(const std::string& k) {
    // common flags across the codebase
    static const std::unordered_set<std::string> whitelist = {
        "print_mpi","numerically_exact","restricted","localize","rotate","checkpoint",
        "dynamic","shared_memory","include_nuclear","include_coulomb","include_xc",
        "isAZORA","spin","nonequilibrium","density","unrestricted","debug","verbose",
        "timings","write_orbitals","plots","plot_density","plot_orbitals","enable",
        "disable","use_gpu","use_mpi","use_omp"
    };

    if (whitelist.count(k)) return true;
    // Heuristics by prefix/suffix
    auto starts_with = [&](const char* p){
        return k.rfind(p, 0) == 0;
    };
    auto ends_with = [&](const char* sfx){
        return k.size() >= std::char_traits<char>::length(sfx)
            && k.compare(k.size()-std::char_traits<char>::length(sfx),
                         std::char_traits<char>::length(sfx), sfx) == 0;
    };

    return starts_with("is_") || starts_with("has_") || starts_with("use_")
        || starts_with("enable_") || starts_with("disable_") || starts_with("do_")
        || starts_with("with_") || starts_with("print_") || starts_with("plot_")
        || ends_with("_enabled") || ends_with("_disabled");
}

// Recursively coerce 0/1 and on/off-like strings to booleans for "boolean-ish" keys.
inline void sanitize_booleans(nlohmann::json& j) {
    if (j.is_object()) {
        for (auto it = j.begin(); it != j.end(); ++it) {
            const std::string key = it.key();
            auto& val = it.value();
            // Recurse first
            if (val.is_object() || val.is_array()) sanitize_booleans(val);

            if (looks_like_bool_key(key)) {
                if (!val.is_boolean()) {
                    // Only coerce if convertible; otherwise leave as-is
                    bool coerced = to_bool(val, /*def*/false);
                    // If val was non-convertible string (e.g., "maybe"), keep old value
                    if (val.is_boolean() || val.is_number() || val.is_string()) {
                        // Special case: avoid mis-coercing numbers not 0/1
                        if (val.is_number_integer()) {
                            int n = val.get<int>();
                            if (n == 0 || n == 1) { val = coerced; }
                        } else if (val.is_number_float()) {
                            double x = val.get<double>();
                            if (x == 0.0 || x == 1.0) { val = coerced; }
                        } else if (val.is_string()) {
                            // Strings like "on"/"off"/"yes"/"no"/"true"/"false"/"0"/"1"
                            std::string s = val.get<std::string>();
                            std::string sl = s;
                            std::transform(sl.begin(), sl.end(), sl.begin(),
                                           [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
                            if (sl=="true"||sl=="false"||sl=="on"||sl=="off"||sl=="yes"||sl=="no"||sl=="0"||sl=="1") {
                                val = coerced;
                            }
                        }
                    }
                }
            }
        }
    } else if (j.is_array()) {
        for (auto& x : j) sanitize_booleans(x);
    }
}

} // namespace mrchem::json_utils
