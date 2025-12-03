from django.shortcuts import render
from .settings import MEDIA_ROOT
from pathlib import Path
import os


DIR_NAME = Path(__file__).resolve().parent


def index(request):
    return render(request, 'index.html')


def analysis(request):
    ddpcr = request.FILES.get('ddpcr')
    ddpcr_file = Path(MEDIA_ROOT / ddpcr.name)
    sheet = request.POST.get('sheet')
    with open(ddpcr_file, 'ab') as fp:
        for chunk in ddpcr.chunks():
            fp.write(chunk)
    run = f'python {DIR_NAME / "ddpcr.py"} -sheet {sheet} -ddpcr {ddpcr_file} -out {ddpcr_file}.analysed.xlsx'
    os.system(run)
    result = {'result': f'{ddpcr.name}.analysed.xlsx'}
    return render(request, 'result.html', result)


def download(request):
    file_name = request.POST.get('file_name')
    result = {'result': file_name}
    return render(request, 'result.html', result)


def transport(request):
    return render(request, 'download.html')
