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
  printf("starting test to build parse and build net form %s", test_case);

  // Run tests
  Parser* parser;
  Net* net;

  printf("test_net_load ... ");

  parser = PARSER_new_for_file(test_case);
  net = PARSER_parse(parser, test_case, 1);

  NET_check(net, 1);
  NET_show_components(net, 2);

  NET_del(net);
  PARSER_del(parser);

  printf("complete");
}
