import time
from ipywidgets import IntProgress, HTML, VBox
from IPython.display import display

# from https://github.com/alexanderkuk/log-progress

def log_progress(sequence, every=None, size=None):
    is_iterator = False
    start_tic = time.time()
    if size is None:
        try:
            size = len(sequence)
        except TypeError:
            is_iterator = True
    if size is not None:
        if every is None:
            if size <= 200:
                every = 1
            else:
                every = size / 200     # every 0.5%
    else:
        assert every is not None, 'sequence is iterator, set every'

    if is_iterator:
        progress = IntProgress(min=0, max=1, value=1)
        progress.bar_style = 'info'
    else:
        progress = IntProgress(min=0, max=size, value=0)
    label = HTML()
    box = VBox(children=[label, progress])
    display(box)

    index = 0
    try:
        for index, record in enumerate(sequence, 1):
            if index == 1 or index % every == 0:
                if is_iterator:
                    label.value = '{index} / ?'.format(index=index)
                else:
                    progress.value = index
                    label.value = u'{index} / {size}'.format(
                        index=index,
                        size=size
                    )
            yield record
    except:
        progress.bar_style = 'danger'
        raise
    else:
        progress.bar_style = 'success'
        progress.value = index
        # pretty human readable time diff
        if False:
            import dateutil.relativedelta
            attrs = ['years', 'months', 'days', 'hours', 'minutes', 'seconds']
            delta = dateutil.relativedelta.relativedelta(seconds=time.time()-start_tic)
            elapsed = " ".join(['%d %s' % (getattr(delta, attr),
                                            getattr(delta, attr) > 1 and attr or attr[:-1]) for attr in attrs if
                                            getattr(delta, attr)])+" %d usec"%((tdiff-int(tdiff))*1000000)
            label.value = u'{index} : {elapsed}'.format(index=index or '?', elapsed=elapsed)
        else:
            # simple time in sec
            label.value = u'{index} : {elapsed:0.2f}s'.format(index=index or '?', elapsed=time.time()-start_tic)

