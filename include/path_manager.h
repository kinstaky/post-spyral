#ifndef __PATH_MANAGER_H__
#define __PATH_MANAGER_H__

#include <mutex>
#include <thread>
#include <filesystem>

namespace gong {

class PathManager {
private:
	PathManager() = default;
	~PathManager() = default;
	PathManager(const PathManager&) = delete;
	PathManager& operator=(const PathManager&) = delete;

public:

	/// @brief Get the Path Manager instance
	/// @returns instance of Path Manager
	static PathManager& Instance();


	/// @brief Initialize the Path Manager, run in main function
	/// @param path argv[0]
	void Initialize(const char *path);


	/// @brief Get the Executable Path
	/// @returns path of executable file
	std::filesystem::path ExecPath() const;


	/// @brief Get the data path
	/// @returns path of data directory
	std::filesystem::path DataPath() const;
private:
	mutable std::mutex mutex_;
	bool init_ = false;
	std::filesystem::path exec_path_;
	std::filesystem::path data_path_;
};


}	// namespace gong


#endif // __PATH_MANAGER_H__