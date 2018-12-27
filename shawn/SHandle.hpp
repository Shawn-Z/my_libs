#ifndef SHAWN_HANDLE_HPP
#define SHAWN_HANDLE_HPP

#include <vector>
#include <ctime>
#include <stdint.h>
#include <stdlib.h>
#include <string>

namespace shawn {
    class handle {
    public:
        friend class SHandle;
    private:
        size_t id;
    };

    class SHandle {
    public:
        handle reserve_handle;

        SHandle() {
            this->id_set.clear();
            this->handle_set.clear();
            this->description_set.clear();
            this->reserve_handle.id = SIZE_MAX;
        }

        handle newHandle(std::string p_handle_description = "none") {
            this->id_set.emplace_back(this->id_set.size());
            handle tmp_handle;
            tmp_handle.id = this->id_set.back();
            this->handle_set.emplace_back(tmp_handle);
            this->description_set.emplace_back(p_handle_description);
            return tmp_handle;
        }

        bool getID(handle p_handle, size_t &p_id) {
            if (p_handle.id > this->id_set.back()) {
                return false;
            }
            p_id = p_handle.id;
            return true;
        }

        bool getDescription(handle p_handle, std::string &p_description) {
            if (p_handle.id > this->id_set.back()) {
                return false;
            }
            p_description = this->description_set[p_handle.id];
            return true;
        }

        inline bool checkHandleCorrect(handle p_handle) {
            return (p_handle.id <= this->id_set.back());
        }

    private:
        std::vector<size_t> id_set;
        std::vector<handle> handle_set;
        std::vector<std::string> description_set;
    };
}

#endif //SHAWN_HANDLE_HPP