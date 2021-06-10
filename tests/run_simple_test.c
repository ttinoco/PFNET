#include <pfnet/parser.h>
#include <pfnet/net.h>

int main(int argc, char **argv) {

  // Local variables
  char* test_case;

  // Check inputs
  if( argc < 2) {
    printf("usage: run_tests test_case\n");
    return -1;
  }

  // Get case
  test_case = argv[1];
  printf("starting test to build parse and build net from %s\n", test_case);

  // Run tests
  Parser* parser;
  Net* net;

  printf("Parsing network ... \n");

  parser = PARSER_new_for_file(test_case);
  PARSER_set(parser, "output_level", 3);
  PARSER_set(parser, "keep_all_out_of_service", FALSE);

  printf("Loading network from parser ... \n");
  net = PARSER_parse(parser, test_case, 1);

  void* opt_out = PARSER_get_setting(parser, "output_level");
  void* opt_oos = PARSER_get_setting(parser, "keep_all_out_of_service");

  if (opt_out)
    printf("output_level = %d\n", *(int*)opt_out);
  if (opt_oos)
    printf("keep_all_out_of_service = %d\n", *(int*)opt_oos);

  NET_check(net, 1);
  NET_show_components(net, 2);

  printf("Cleaning up ... \n");

  NET_del(net);
  PARSER_del(parser);

  printf("complete!\n");
}
