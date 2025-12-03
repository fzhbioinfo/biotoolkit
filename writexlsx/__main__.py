"""解读表格汇总
"""
import json
from pathlib import Path
from argparse import ArgumentParser
from writexlsx.utils import WriteXlsx


BASE_DIR = Path(__file__).resolve().parent


def run(args):
    """解读表格汇总

    :param args: Namespace
    :type args: Namespace
    """
    with open(args.json, 'r', encoding='utf-8') as j:
        result_dic = json.load(j)
    template = args.template if args.template else BASE_DIR / f'template/template_{result_dic["template"]}.xlsx'
    with open(args.config, 'r', encoding='utf-8') as j:
        conf_dic = json.load(j)
    try:
        style = conf_dic[result_dic["template"]]
    except KeyError:
        style = {}
    write_xlsx = WriteXlsx(template, style)
    for sheet_name in result_dic:
        if sheet_name == 'template':
            continue
        if result_dic[sheet_name].endswith('.json'):
            write_xlsx.write_json_to_sheet(sheet_name, result_dic[sheet_name])
        elif result_dic[sheet_name].endswith('.list'):
            write_xlsx.write_list_to_sheet(sheet_name, result_dic[sheet_name])
        else:
            write_xlsx.write_table_to_sheet(sheet_name, result_dic[sheet_name])
    write_xlsx.save_book(args.out)


def main():
    """模块输入参数
    """
    parser = ArgumentParser(description='Write variant detected results into excel')
    parser.add_argument('-j', dest='json', help='result json', required=True)
    parser.add_argument('-o', dest='out', help='output file name', required=True)
    parser.add_argument('-t', dest='template', help='excel template')
    parser.add_argument('-c', dest='config', default=BASE_DIR / 'config/style.json',
                        help='style config')
    parsed_args = parser.parse_args()
    run(parsed_args)


if __name__ == '__main__':
    main()
