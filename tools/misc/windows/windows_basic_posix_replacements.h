// ----------------------------------------------------------------------------------------
// 	The purpose of this file: Osiris uses a few functionswhich are only found in POSIX. 
//	  We need to provide funcitons that map the POSIX functions to the corresponding Windows ones.
// ----------------------------------------------------------------------------------------

#ifndef _WINDOWS_POSIX_REPLACEMENTS_
#define _WINDOWS_POSIX_REPLACEMENTS_
	
	
	#define _CRT_SECURE_NO_WARNINGS 1
	#pragma warning(disable:4996)
	//
	// Include needed headers for replacement functions. Serveral are Windows-only.
	//
	//		io.h - needed for the 'umask' functions.
	#include <windows.h>
	#include <direct.h>
	#include <float.h>
	#include <tchar.h>
	#include <time.h>
	#include <io.h>
	#include <stdlib.h>
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <errno.h>
	#include <string.h>
	#include <stdio.h>
	#include <math.h>

	//
	// Set Windows API configuration variables used in this file.
	//
	// 		WIN_HOST_NAME_TYPE = 3 is a Windows enumeration to report local NETBIOS host name
	#define WIN_HOST_NAME_TYPE 3	
	

	// ----------------------------------------------------------------------------------------
	//	Several of the missing POSIX functions in Windows can
	//	  be replaced by simple renaming to their Windows cognate.
	// ----------------------------------------------------------------------------------------
	
	//	The windows "mkdir" is actually more standards compliant and take only 1 argument
	//		So we simply ignore the 2nd argument passed in.
	#define mkdir(a,b) _mkdir(a)
	
	// it's not necessary to do these mappings.. but it gets rid of som visual studio warnings.
	//		what's happening is that the Microsoft compiler is actually more vocal about standards
	//		compliance .. and technically chdir, mkdir, getcwd are depricated in favor of
	//		_chdir, _mkdir, and _getcwd
	#define chdir(a) _chdir(a)
	#define getcwd _getcwd
	#define umask(a) _umask(a)
	
	
	// ----------------------------------------------------------------------------------------
	// Windows has some wonky, non standard math stuff.. some of it is better alot is not
	// ----------------------------------------------------------------------------------------
		
	//
	// Visual Studio 2008 and earlier need replacement definitions for some basic math functions.
	// Visual studio 2013 standardized most of the math functions (adding most of C99) 
	// and the replacements are not needed.
	//
	// I tell which case we are in by using the fact that 2013 
	// defined 'isnan' and  'isinf' as macros so I can look for those. 
	// It's a bit of fact.. should select on _MSC_VER variable.
	// TODO: do this the official nice way.
	//
	// For VS 2008 and earlier, the replacement functions are:
	#if !defined( isnan ) &&  !defined( isinf )
		
		#define isnan(x) _isnan(x)
		#define isinf(x) (!_finite(x) && !_isnan(x))
		
		#define INFINITY_D (DBL_MAX+DBL_MAX)
		#define NAN_D (INFINITY_D-INFINITY_D)
		#define INFINITY_F (FLT_MAX+FLT_MAX)
		#define NAN_F (INFINITY_F-INFINITY_F)

		double nan( const char * x) { return (NAN_D); }
		float nanf( const char * x) { return (NAN_F); }
	#endif	
	
	// Just a quick test function for the replaced floating point stuff
	void quick_windows_FP_test() {
		double NAND; double NANF;
		NAND = nan("");
		printf("isnan(NAND) %d \n", isnan(NAND) );
		printf("isinf(NAND) %d \n", isinf(NAND) );
		printf("isfinite(NAND) %d \n", _finite(NAND) );
		NANF = nanf("");
		printf("isnan(NANF) %d \n", isnan(NANF) );
		printf("isinf(NANF) %d \n", isinf(NANF) );
		printf("isfinite(NANF) %d \n", _finite(NANF) );
	}

	
	// ----------------------------------------------------------------------------------------
	//	Now drop in a replacement for the POSIX timing functions...
	// ----------------------------------------------------------------------------------------

	// a simple forward declaration.. it needed in Windows for various reasons
	double MPI_Wtime( void );
	double MPI_Wtick( void );
		
	#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
		#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
	#else
		#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
	#endif

	struct timezone
	{
	  int  tz_minuteswest; /* minutes W of Greenwich */
	  int  tz_dsttime;     /* type of dst correction */
	};
	
	
	int gettimeofday(struct timeval *tv, struct timezone *tz) {
	  
	  // Define a structure to receive the current Windows filetime
	  FILETIME ft;
	 
		// Initialize the present time to 0 and the timezone to UTC
	  unsigned __int64 tmpres = 0;
	  static int tzflag = 0;
	 
	  if (NULL != tv) {
		GetSystemTimeAsFileTime(&ft);
	 
		// The GetSystemTimeAsFileTime returns the number of 100 nanosecond 
		// intervals since Jan 1, 1601 in a structure. Copy the high bits to 
		// the 64 bit tmpres, shift it left by 32 then or in the low 32 bits.
		tmpres |= ft.dwHighDateTime;
		tmpres <<= 32;
		tmpres |= ft.dwLowDateTime;
	 
		// Convert to microseconds by dividing by 10
		tmpres /= 10;
	 
		// The Unix epoch starts on Jan 1 1970.  Need to subtract the difference 
		// in seconds from Jan 1 1601.
		tmpres -= DELTA_EPOCH_IN_MICROSECS;
	 
		// Finally change microseconds to seconds and place in the seconds value. 
		// The modulus picks up the microseconds.
		tv->tv_sec = (long)(tmpres / 1000000UL);
		tv->tv_usec = (long)(tmpres % 1000000UL);
	  }
	 
	  if (NULL != tz) {
		if(!tzflag) {
			_tzset();
			tzflag++;
		}
		// Adjust for the timezone west of Greenwich
		//tz->tz_minuteswest = _timezone / 60;
		//tz->tz_dsttime = _daylight;
	  }
	  return 0;
	}
	
	// ----------------------------------------------------------------------------------------
	//	Replacement code for creating Symbolic Links
	// ----------------------------------------------------------------------------------------
	int symlink(const char *name1, const char *name2) {
		// check if are in Unicode mode.. then we have some conversion to do
		//		cause in Unicode mode, the  Windows API function that symbolic links 
		//		wants unicode strings. So the simple char* fortran input must be 
		//		converted to unicode ( using UTF-8 encoding )
		#if defined(UNICODE) || defined(_UNICODE)
			int size_needed;
			// WCHAR is always 16-bit char,
			// TCHAR changes from CHAR to WCHAR depending on UNICODE #defines.. 
			//		TCHAR is the cooler way to do things..
			// 		But we just use WCHAR since this code is only 
			//		compiled if Unicode is enabled.
			WCHAR unicode_buffer__link[1024];	
			WCHAR unicode_buffer__target[1024];

			size_needed = MultiByteToWideChar(CP_UTF8, 0, &name1[0], -1, NULL, 0);
			if( size_needed > 1024) {
				printf("Error creating symbolic link: link path too long: %s\n", name1);
				return -1;
			}
			MultiByteToWideChar(CP_UTF8, 0, &name1[0], -1, &unicode_buffer__link[0], size_needed);
			
			size_needed = MultiByteToWideChar(CP_UTF8, 0, &name2[0], -1, NULL, 0);
			if( size_needed > 1024) {
				printf("Error creating symbolic link: target path too long: %s\n", name2);
				return -1;
			}
			MultiByteToWideChar(CP_UTF8, 0, &name2[0], -1, &unicode_buffer__target[0], size_needed);
			
			// LPCTSTR is a LPCWSTR (const 16 bit char pointer) if UNICODE defined.. 
			//		otherwise it's a LPCSTR (const 8 bit char pointer )
			LPCTSTR link = name2;
			LPCTSTR target = name1;
		#else
			// This is the usual case.. we are NOT in Unicode mode on the windows API side
			//		and the incoming char* is good as-is.
			LPCTSTR link = name2;
			LPCTSTR target = name1;			
		#endif
		#undef UNICODE
		
		BOOL fCreatedLink;
		// Works only in vista and above (which is OK these days)
		//
		// note: the following works in all windows versions:
		// 			fCreatedLink = CreateHardLink(link,realFile,NULL);
		// 		But on there is a hitch on Windows XP.. 
		//		XP can't handle relative paths in symbolic links.....
		// 		If full XP support is desired, then:
		//			find way to get current working directory and prepend it
		//			if XP support ever becomes necessary.
		//		0x2 is 'SYMBOLIC_LINK_FLAG_ALLOW_UNPRIVILEGED_CREATE'. For Windows 10 and above, it allows
		//			anyuser to make sym links (in earlier windows versions, you had to run Osiris as Administrator 
		//			for the restart symlinks to work)
		fCreatedLink = CreateSymbolicLink(link, target, (DWORD) 0x2 );

		if(fCreatedLink == FALSE ) {
			// error creating the symbolic link...
			return errno;
		}
		
		// successful return..
		return 0;
	}
	
	// ----------------------------------------------------------------------------------------
	//	Replacement code for getting Hostname
	//	( No longer needed for VS 2008 and greater).
	// ----------------------------------------------------------------------------------------
	void gethostname_windows( char* hostname, int num_num_characters) {
		// The windows way to get the hostname
		// windows magic type buffer to hold the buffer..
		TCHAR buffer[256];
		DWORD dwSize = sizeof(buffer);
		
		if(!GetComputerNameEx((COMPUTER_NAME_FORMAT) WIN_HOST_NAME_TYPE, buffer, &dwSize)) {
			// host name lookup failed..
			hostname[0] = 'u';
			hostname[1] = 'n';
			hostname[2] = 'k';
			hostname[3] = 'o';
			hostname[4] = 'w';
			hostname[5] = 'n';
			hostname[6] = '\0';
		} else {
			strcpy(hostname,buffer);
		}
		
		// ensure string is terminated.
		strcpy(&hostname[255],"\0");
	}

  void does_dir_exist_f(const char * path, int* exists)
  {
    DWORD attribs = GetFileAttributesA(path);
    if (attribs == INVALID_FILE_ATTRIBUTES) {
      /* Directory does not exist */
      *exists = 0;
      
    } else {
      *exists = 1;
    }
  }

  /*
    Get the number of bytes in file 'filename'
      returns -1 if an error occured.
      returns -1 if the file size exceeds 2^31-1 bytes (i.e. ~ 2 gigabytes)
  */
  int32_t get_filesize(const TCHAR *filename) {
    #define MAX_CRC_BUFFER_LEN 262144
    BOOL                        is_ok;
    WIN32_FILE_ATTRIBUTE_DATA   file_info;
    is_ok = GetFileAttributesEx(filename, GetFileExInfoStandard, (void*)&file_info);
    if (!is_ok) return -1;
    
    if(file_info.nFileSizeHigh != 0) return -1;
    if(file_info.nFileSizeLow >= INT32_MAX) return -1;
    return (uint32_t) file_info.nFileSizeLow;
  }

#endif