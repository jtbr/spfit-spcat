#include <time.h>
#include <string.h>

int strdate_(char *pbuffer, int plen)
{
  time_t curtime;
  struct tm *loctime;
  char temp_buffer[64]; /* Temporary buffer for strftime output */
  size_t result_len;
  int i;

  /* Get the current time. */
  time(&curtime);

  /* Convert it to local time representation. */
  loctime = localtime(&curtime);

  /* Format the date and time using strftime - format matches asctime output */
  result_len = strftime(temp_buffer, sizeof(temp_buffer), "%a %b %d %H:%M:%S %Y", loctime);

  /* Copy to output buffer with padding if needed */
  for (i = 0; i < plen; i++) {
    if (i < result_len) {
      pbuffer[i] = temp_buffer[i];
    } else {
      pbuffer[i] = ' ';
    }
  }

  return 0;
}
