#ifndef _prapare_directory_hpp_
#define _prapare_directory_hpp_

#include <experimental/filesystem>
#include <stdexcept>

// prepare a directory given by an input string and return the directory_entry object
inline std::experimental::filesystem::directory_entry prepare_directory(std::string path, const bool keep_old_files=true)
{
  std::cout << "\n..." << std::endl;

  if (path[0]=='~'){
    std::string home_dir = getenv("HOME");
    path.erase(0,1);
    path = home_dir + path;
  }

  namespace fs = std::experimental::filesystem;
  fs::directory_entry directory;

  directory.assign(path);
  std::cout << "preparing directory: " << directory.path() << std::endl;

  if (not fs::exists(directory.path()))
  {
    // std::cout << "warning: output directory does NOT exist!!!" << std::endl;
    std::cout << "created the directory!" << std::endl;
    fs::create_directories(directory.path());
    if (not fs::is_directory(directory.path())) throw std::invalid_argument("The input value for output directory is not acceptable.");
    return directory;
  }

  if (fs::is_directory(directory.path()))
  {
    if (not fs::is_empty(directory.path()))
    {
      std::cout << "warning: output directory is NOT empty!!!" << std::endl;
      if (keep_old_files)
      {
        int count = 1;
        while (fs::exists(directory.path().string()+"."+std::to_string(count))){
          count++;
        }
        fs::path new_path = directory.path().string()+"."+std::to_string(count);
        std::cout << "renaming the existing directory to " << new_path << std::endl;
        fs::rename(directory.path(),new_path);
      } 
      else
      {
        std::cout << "deleting the existing directory!!!" << std::endl;
        fs::remove_all(directory.path());
      }
      fs::create_directories(directory.path());
    }
  }
  else
  {
    throw std::invalid_argument("The input value for output directory is not acceptable.");
  }

  std::cout << "...\n" << std::endl;

  return directory;

};

// check a string input to make sure it points to an existing directory and return the directory_entry
inline std::experimental::filesystem::directory_entry check_directory(std::string path, bool should_be_empty=false)
{

  std::cout << "\n..." << std::endl;

  if (path[0]=='~'){
    std::string home_dir = getenv("HOME");
    path.erase(0,1);
    path = home_dir + path;
  }

  namespace fs = std::experimental::filesystem;
  fs::directory_entry directory;

  directory.assign(path);
  std::cout << "checking directory: " << directory.path() << std::endl;

  if (not fs::exists(directory.path()))
  {
    throw std::invalid_argument("directory does NOT exists!!!");
  }

  if (fs::is_directory(directory.path()))
  {
    if (fs::is_empty(directory.path()) and (not should_be_empty))
    {
      throw std::invalid_argument("directory is empty!!!");
    }
  }
  else
  {
    throw std::invalid_argument("input path is NOT a directory!!!");
  }
  
  std::cout << "...\n" << std::endl;

  return directory;
}

#endif //_prapare_directory_hpp_