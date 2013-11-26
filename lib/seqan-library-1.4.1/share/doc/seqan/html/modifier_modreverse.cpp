#include <iostream>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace seqan;


int main ()
{
	String<char> myString = "A man, a plan, a canal-Panama";
	ModifiedString< String<char>, ModReverse > myModifier(myString);

	std::cout << myString << std::endl;
	std::cout << myModifier << std::endl;
	replace(myString, 9, 9, "master ");
	std::cout << myString << std::endl;
	std::cout << myModifier << std::endl;
	return 0;
}
